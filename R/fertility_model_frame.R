#' Aggregate district level populations
#' @description Aggregate district level populations to all levels of area hierarchy
#' @param population Age/sex/space/time stratified population
#' @param areas_wide Area hierarchy in wide format
#' @export

area_populations <- function(population, areas_wide, max_year) {

  population <- population %>%
    filter(sex == "female", period <= max_year)

  base_area_pop <- left_join(areas_wide %>% st_drop_geometry(), population, by = c("area_id"))

  level_ids <- str_subset(colnames(areas_wide), "area_id[0-9]")

  group_area_pops <- function(level_ids, base_area_pop) {

    base_area_pop %>%
      group_by(.data[[level_ids]], sex, age_group, period) %>%
      summarise(population = sum(population),
                .groups = "drop") %>%
      rename(area_id = .data[[level_ids]])

  }

  area_populations <- lapply(level_ids, group_area_pops, base_area_pop) %>%
    bind_rows %>%
    filter(!is.na(area_id))

  return(area_populations)

}

#' Make random walk structure matrices
#' @description Make random walk structure matrices
#' @param x Matrix size
#' @param order Random walk order
#' @param adjust_diagonal Add 1E-6 to the matrix diagonal to make the matrix proper. Default = TRUE

make_rw_structure_matrix <- function(x, order, adjust_diagonal = TRUE) {



  D_mat <- diff(diag(x), differences = order)
  R_mat <- t(D_mat) %*% D_mat

  if(adjust_diagonal) {
    diag(R_mat) <- diag(R_mat) + 1E-6
  }

  R_mat <- as(R_mat, "dgCMatrix")

  return(R_mat)

}

make_adjacency_matrix <- function(areas) {

  sh <- areas %>%
    filter(naomi_level) %>%
    arrange(area_sort_order) %>%
    mutate(area_idx = row_number())

  #' Neighbor list
  nb <- sh %>%
    as("Spatial") %>%
    spdep::poly2nb() %>%
    `names<-`(sh$area_idx)

  names(nb) <- sh$area_id


  nb <- lapply(nb, as.integer)
  class(nb) <- "nb"

  adj <- spdep::nb2mat(nb, zero.policy=TRUE, style="B")
  R_spatial <- INLA::inla.scale.model(diag(rowSums(adj)) - 0.99*adj,
                                      constr = list(A = matrix(1, 1, nrow(adj)), e = 0))

  return(R_spatial)

}

#' Make model frames
#' @description Create model frames and aggregation matrices for TMB model.
#' @iso3
#' @param population Age/sex/space stratified population
#' @param asfr ASFRs by district and time
#' @param mics_asfr ASFRs at provincial level from MICS surveys
#' @param areas Area hierarchy
#' @param project Model will by default produce estimates up to year of last survey. Integer projection year if desired, else FALSE
#' @export
#'
make_model_frames <- function(iso3_c,
                              population,
                              asfr,
                              mics_asfr,
                              areas,
                              project = FALSE) {

  areas_long <- areas %>%
    st_drop_geometry() %>%
    mutate(iso3 = iso3_c)

  if(is.null(mics_asfr))
    mics_flag <- FALSE
  else
    mics_flag <- TRUE

  asfr <- asfr %>%
    bind_rows()

  if(!project) {
    if(mics_flag) {
      mics_asfr <- mics_asfr %>%
        bind_rows()

      df <- asfr %>%
        bind_rows(mics_asfr)
    } else {
      df <- asfr
    }

    max_year <- max(df$period)

  } else {
    max_year <- project
  }

  population <- area_populations(population, spread_areas(areas), max_year) %>%
    filter(sex == "female") %>%
    ungroup %>%
    dplyr::select(-sex)

  population <- crossing(area_id = unique(population$area_id),
                         age_group = unique(population$age_group),
                         period = 1995:2020
  ) %>%
    left_join(population, c("area_id", "age_group", "period")) %>%
    group_by(area_id, age_group) %>%
    mutate(population = exp(zoo::na.approx(log(population), period, na.rm = FALSE))) %>%
    ungroup() %>%
    fill(population, .direction="updown") %>%
    left_join(
      areas_long %>% dplyr::select(area_id, iso3),
      by = "area_id"
    )

  area_tree <- create_areas(area_merged = areas)
  area_aggregation <- create_area_aggregation(areas$area_id[areas$naomi_level], area_tree)

  mf_model <- crossing(period = 1995:max_year,
                       age_group = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
                       area_id = unique(area_aggregation$model_area_id)) %>%
    left_join(
      population %>%
        dplyr::select(iso3, area_id, period, age_group, population),
      by = c("period", "age_group", "area_id")
    ) %>%
    mutate(area_id = factor(area_id),
           age_group = factor(age_group, levels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
           period = factor(period),
           urban = ifelse(area_id %in% c(
             filter(areas_long, parent_area_id == "ETH_1_10")$area_id,
             filter(areas_long, str_detect(area_name, "Town"))$area_id,
             filter(areas_long, area_name %in% c("Harari", "Fafen (Jijiga)"))$area_id),
             1, 0)
    ) %>%
    arrange(period, area_id, age_group) %>%
    mutate(idx = factor(row_number()),
           id.interaction1 = factor(group_indices(., age_group, period, iso3)),
           id.interaction2 = factor(group_indices(., period, area_id)),
           id.interaction3 = factor(group_indices(., age_group, area_id)),
           id.omega1 = factor(group_indices(., age_group, iso3)),
           id.omega2 = factor(group_indices(., period, iso3))
    )

  obs <- asfr %>%
    mutate(period = factor(period, levels(mf_model$period))) %>%
    filter(!is.na(births)) %>%
    dplyr::select(survtype, area_id, period, age_group, tips, births, pys) %>%
    left_join(mf_model, by = c("area_id", "period", "age_group")) %>%
    mutate(tips_dummy = as.integer(tips > 5),
           tips_f = factor(tips),
           ais_dummy = ifelse(survtype %in% c("MIS", "AIS"), 1, 0),
           #####
           urban = ifelse(area_id %in% c(
             filter(areas_long, parent_area_id == "ETH_1_10")$area_id,
             filter(areas_long, str_detect(area_name, "Town"))$area_id,
             filter(areas_long, area_name %in% c("Harari", "Fafen (Jijiga)"))$area_id),
             1, 0),
           #####
           age_group = factor(age_group, levels(mf_model$age_group)),
           area_id = factor(area_id, levels(mf_model$area_id)),
           period = factor(period, levels(mf_model$period)),
    )

  mf <- list()
  mf$mf_model <- mf_model
  mf$district$obs <- obs
  mf$mics_toggle <- 0
  mf$out_toggle <- 0

  ## Outputs

  age_aggregation <- data.frame(
    "age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group_label,
    "model_age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group_label
  ) %>%
    bind_rows(data.frame(
      "age_group" = "15-49",
      "model_age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group_label
    ))

  asfr_out <- crossing(
    area_id = area_aggregation$area_id,
    age_group = unique(mf_model$age_group),
    period = unique(mf_model$period),
    variable = "asfr"
  ) %>%
    arrange(variable, area_id, age_group, period) %>%
    mutate(out_idx = row_number()) %>%
    droplevels()

  asfr_join_out <- crossing(area_aggregation,
                            age_group = unique(mf_model$age_group),
                            period = unique(mf_model$period)) %>%
    full_join(mf_model %>%
                dplyr::select(area_id, age_group, period, idx),
              by = c("model_area_id" = "area_id", "age_group", "period")
    ) %>%
    full_join(asfr_out, by = c("area_id", "age_group", "period")) %>%
    mutate(x=1) %>%
    filter(!is.na(model_area_id))

  tfr_out <- crossing(
    area_id = area_aggregation$area_id,
    age_group = "15-49",
    period = unique(mf_model$period),
    variable = "tfr"
  ) %>%
    arrange(variable, area_id, age_group, period) %>%
    mutate(out_idx = row_number()) %>%
    droplevels()

  asfr_join_out %>%
    dplyr::select(area_id, age_group, period, out_idx) %>%
    unique()

  tfr_join_out <- crossing(area_id = area_aggregation$area_id,
                           age_aggregation %>% filter(age_group == "15-49"),
                           period = unique(mf_model$period)) %>%
    full_join(
      asfr_join_out %>%
        dplyr::select(area_id, age_group, period, out_idx) %>%
        unique,
      by = c("area_id", "model_age_group" = "age_group", "period")
    ) %>%
    rename(idx = out_idx) %>%
    arrange(period, area_id, age_group) %>%
    full_join(tfr_out, by = c("area_id", "age_group", "period")) %>%
    mutate(x=5)

  A_asfr_out <- spMatrix(nrow(asfr_out), nrow(mf_model), asfr_join_out$out_idx, as.integer(asfr_join_out$idx), asfr_join_out$x)
  A_tfr_out <- spMatrix(nrow(tfr_out), nrow(asfr_out), tfr_join_out$out_idx, as.integer(tfr_join_out$idx), tfr_join_out$x)

  mf_out <- asfr_out %>%
    bind_rows(tfr_out) %>%
    dplyr::select(-out_idx)

  mf$out$mf_out <- mf_out
  mf$out$A_asfr_out <- A_asfr_out
  mf$out$A_tfr_out <- A_tfr_out
  mf$out_toggle <- 1


  if(mics_flag) {

    # mf_mics <- crossing(area_id = filter(areas_long, area_level == 1)$area_id,
    mf_mics <- crossing(area_id = as.character(unique(mics_asfr$area_id)),
                        period = unique(mf_model$period),
                        age_group = unique(mf_model$age_group)
    ) %>%
      mutate(idx = factor(row_number()))

    join_mics <- mf_mics %>%
      rename(idx_row = idx) %>%
      left_join(area_aggregation, by = "area_id") %>%
      left_join(mf_model, by=c("age_group", "period", "model_area_id" = "area_id")) %>%
      mutate(x=1) %>%
      type.convert()

    A_mics <- sparseMatrix(i = join_mics$idx_row, j=join_mics$idx, x=join_mics$x, dims = c(max(join_mics$idx_row), nrow(mf_model)), use.last.ij = TRUE)


    obs_mics <- mics_asfr %>%
      mutate(period = factor(period, levels(mf_model$period))) %>%
      left_join(mf_mics, by = c("area_id", "age_group", "period")) %>%
      dplyr::select(area_id, period, age_group, tips, births, pys, idx) %>%
      mutate(tips_dummy = as.integer(tips > 2),
             tips_f = factor(tips),
             age_group = factor(age_group, levels(mf_model$age_group)),
             idx =factor(idx, levels(mf_mics$idx))
      )

    mf$mics$obs <- obs_mics
    mf$mics$A_mics <- A_mics
    mf$mics_toggle <- 1

  }

  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  Z$Z_country <- as(matrix(rep(1, nrow(mf$mf_model)), ncol=1), "dgCMatrix")

  M_obs <- sparse.model.matrix(~0 + idx, mf$district$obs)
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$district$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$district$obs %>% filter(ais_dummy ==0))
  # Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$district$obs %>% filter(ais_dummy ==1))
  # Z_tips[which(mf$district$obs$survtype != "DHS"), ] <- 0
  Z$X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$district$obs %>% filter(ais_dummy == 0))
  Z$X_urban_dummy <- model.matrix(~0 + urban, mf$district$obs)

  ais_join <- mf$district$obs %>%
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)

  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$district$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)

  dhs_join <- mf$district$obs %>%
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)

  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$district$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)

  X_extract <- list()
  X_extract$X_extract_ais <- X_extract_ais
  X_extract$X_extract_dhs <- X_extract_dhs

  R <- list()
  R$R_spatial <- make_adjacency_matrix(areas)
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  R$R_country <- as(diag(1, 1), "dgTMatrix")
  R$R_spatial_iid <- as(diag(1, length(unique(mf$mf_model$area_id))), "dgTMatrix")

  if(mics_flag) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs)
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    # X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
    R$R_tips_mics <- make_rw_structure_matrix(ncol(Z$Z_tips_mics), 1, adjust_diagonal = TRUE)
  }

  mf$R <- R
  mf$Z <- Z
  mf$X_extract <- X_extract

  mf$district$M_obs <- M_obs
  mf$mics$M_obs_mics <- M_obs_mics

  return(mf)
}
