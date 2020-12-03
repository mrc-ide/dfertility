#' Aggregate district level populations
#' @description Aggregate district level populations to all levels of area hierarchy
#' @param population Age/sex/space/time stratified population
#' @param areas_wide Area hierarchy in wide format
#' @export

area_populations <- function(population, areas_wide, project) {

  population <- population %>%
    dplyr::filter(sex == "female", period <= project)

  base_area_pop <- dplyr::left_join(
    areas_wide %>% sf::st_drop_geometry(), population, by = c("area_id"))

  level_ids <- stringr::str_subset(colnames(areas_wide), "area_id[0-9]")

  group_area_pops <- function(level_ids, base_area_pop) {

    base_area_pop %>%
      dplyr::group_by(.data[[level_ids]], sex, age_group, period) %>%
      dplyr::summarise(population = sum(population),
                .groups = "drop") %>%
      dplyr::rename(area_id = .data[[level_ids]])

  }

  area_populations <- lapply(level_ids, group_area_pops, base_area_pop) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!is.na(area_id))

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

make_adjacency_matrix <- function(areas, model_level) {

  sh <- areas %>%
    filter(area_level == model_level) %>%
    arrange(area_sort_order) %>%
    mutate(area_idx = row_number())

  #' Neighbor list
  nb <- sh %>%
    # st_as_sf() %>%
    # as("Spatial") %>%
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


make_model_frames_dev <- function(iso3,
                              population,
                              asfr,
                              areas,
                              naomi_level,
                              project = 2020) {

  areas_long <- areas %>%
    st_drop_geometry() %>%
    mutate(iso3 = iso3)

  population <- population %>%
    mutate(period = as.numeric(year_labels(calendar_quarter_to_quarter_id(calendar_quarter)))) %>%
    select(-calendar_quarter)

  population <- area_populations(population, spread_areas(areas), project) %>%
    filter(sex == "female") %>%
    ungroup %>%
    dplyr::select(-sex)

  population <- tidyr::crossing(area_id = unique(population$area_id),
                         age_group = unique(population$age_group),
                         period = 1995:2020
  ) %>%
    left_join(population, c("area_id", "age_group", "period")) %>%
    group_by(area_id, age_group) %>%
    mutate(population = exp(zoo::na.approx(log(population), period, na.rm = FALSE))) %>%
    tidyr::fill(population, .direction="updown") %>%
    left_join(
      areas_long %>% dplyr::select(area_id, iso3),
      by = "area_id"
    )

  area_tree <- create_areas(area_merged = areas)
  area_aggregation <- create_area_aggregation(areas$area_id[areas$area_level == naomi_level], area_tree)

  mf_model <- tidyr::crossing(period = 1995:project,
                       age_group = unique(asfr$age_group),
                       area_id = unique(area_aggregation$model_area_id)) %>%
    left_join(
      population %>%
        dplyr::select(iso3, area_id, period, age_group, population),
      by = c("period", "age_group", "area_id")
    ) %>%
    mutate(area_id = factor(area_id),
           age_group = factor(age_group, levels = unique(asfr$age_group)),
           period = factor(period),
           urban = ifelse(area_id %in% c(
             filter(areas_long, parent_area_id == "ETH_1_10")$area_id,
             filter(areas_long, str_detect(area_name, "Town"))$area_id,
             filter(areas_long, area_name %in% c("Harari", "Fafen (Jijiga)", "Erer"))$area_id),
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
    left_join(areas_long) %>%
    mutate(period = factor(period, levels(mf_model$period))) %>%
    filter(!is.na(births)) %>%
    dplyr::select(survtype, area_id, area_level, period, age_group, tips, births, pys) %>%
    # left_join(full_frame) %>%
    mutate(tips_dummy = as.integer(tips > 5),
           tips_f = factor(tips),
           ais_dummy = ifelse(survtype %in% c("MIS", "AIS"), 1, 0),
           mics_dummy = ifelse(survtype == "MICS", 1, 0),
           age_group = factor(age_group, levels(mf_model$age_group)),
           period = factor(period, levels(mf_model$period)),
           spike_2000 = ifelse(period == 2000, 1, 0),
           spike_1999 = ifelse(period == 1999, 1, 0),
           spike_2001 = ifelse(period == 2001, 1, 0)
    )


  naomi_level_obs <- obs %>%
    filter(area_level == naomi_level) %>%
    mutate(
           area_id = factor(area_id, levels(mf_model$area_id))
    ) %>%
    left_join(mf_model %>% select(area_id, age_group, period, idx))

  aggregate_mf <- tidyr::crossing(area_id = area_aggregation$area_id,
                         period = unique(mf_model$period),
                         age_group = unique(mf_model$age_group)
  ) %>%
    arrange(area_id, age_group, period) %>%
    mutate(idx_row = factor(row_number()))

  join <- aggregate_mf %>%
    left_join(area_aggregation, by = "area_id") %>%
    left_join(mf_model, by=c("age_group", "period", "model_area_id" = "area_id")) %>%
    mutate(x=1) %>%
    type.convert()

  full_obs <- obs %>%
    left_join(aggregate_mf, by=c("area_id", "age_group", "period")) %>%
    rename(idx = idx_row)

  A_full_obs <- sparseMatrix(i = join$idx_row, j=join$idx, x=join$x, dims = c(max(join$idx_row), nrow(mf_model)), use.last.ij = TRUE)

  mf <- list()
  mf$mf_model <- mf_model
  mf$observations$naomi_level_obs <- naomi_level_obs
  mf$observations$full_obs <- full_obs
  mf$observations$A_full_obs <- A_full_obs

  ## Outputs

  age_aggregation <- data.frame(
    "age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group,
    "model_age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
  ) %>%
    bind_rows(data.frame(
      "age_group" = "Y015-Y049",
      "model_age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
    ))

  # asfr_out <- tidyr::crossing(
  #   area_id = area_aggregation$area_id,
  #   age_group = unique(mf_model$age_group),
  #   period = unique(mf_model$period),
  #   variable = "asfr"
  # ) %>%
  #   arrange(variable, area_id, age_group, period) %>%
  #   mutate(out_idx = row_number()) %>%
  #   droplevels()
  #
  # asfr_join_out <- tidyr::crossing(area_aggregation,
  #                           age_group = unique(mf_model$age_group),
  #                           period = unique(mf_model$period)) %>%
  #   full_join(mf_model %>%
  #               dplyr::select(area_id, age_group, period, idx),
  #             by = c("model_area_id" = "area_id", "age_group", "period")
  #   ) %>%
  #   full_join(asfr_out, by = c("area_id", "age_group", "period")) %>%
  #   mutate(x=1) %>%
  #   filter(!is.na(model_area_id))

  tfr_out <- crossing(
    area_id = area_aggregation$area_id,
    age_group = "Y015-Y049",
    period = unique(mf_model$period),
    variable = "tfr"
  ) %>%
    arrange(variable, area_id, age_group, period) %>%
    mutate(idx_out = row_number()) %>%
    droplevels()

  tfr_join_out <- crossing(area_id = area_aggregation$area_id,
                           age_aggregation %>% filter(age_group == "Y015-Y049"),
                           period = unique(mf_model$period)) %>%
    full_join(
      aggregate_mf %>%
        dplyr::select(area_id, age_group, period, idx_row) %>%
        unique,
      by = c("area_id", "model_age_group" = "age_group", "period")
    ) %>%
    rename(idx_col = idx_row) %>%
    arrange(period, area_id, age_group) %>%
    full_join(tfr_out, by = c("area_id", "age_group", "period")) %>%
    mutate(x=5)

  # A_asfr_out <- spMatrix(nrow(aggregate_mf), nrow(mf_model), asfr_join_out$out_idx, as.integer(asfr_join_out$idx), asfr_join_out$x)
  A_tfr_out <- spMatrix(nrow(tfr_out), nrow(aggregate_mf), tfr_join_out$idx_out, as.integer(tfr_join_out$idx_col), tfr_join_out$x)

  mf_out <- aggregate_mf %>%
    mutate(variable = "asfr") %>%
    bind_rows(tfr_out) %>%
    dplyr::select(-idx_out)

  mf$out$mf_out <- mf_out
  # mf$out$A_asfr_out <- A_asfr_out
  mf$out$A_tfr_out <- A_tfr_out

  M_naomi_obs <- sparse.model.matrix(~0 + idx, mf$observations$naomi_level_obs)
  M_full_obs <- sparse.model.matrix(~0 + idx, mf$observations$full_obs)

  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  Z$Z_country <- as(matrix(rep(1, nrow(mf$mf_model)), ncol=1), "dgCMatrix")

  # Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$observations$full_obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% filter(survtype == "DHS"))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% filter(survtype %in% c("AIS", "MIS")))
  Z$X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$observations$full_obs %>% filter(survtype == "DHS"))
  Z$X_urban_dummy <- model.matrix(~0 + urban, mf$mf_model)

  ais_join <- mf$observations$full_obs %>%
    mutate(col_idx = row_number()) %>%
    select(col_idx, survtype) %>%
    filter(survtype %in% c("AIS", "MIS")) %>%
    mutate(row_idx = row_number(),
           x=1)

  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$observations$full_obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)

  dhs_join <- mf$observations$full_obs %>%
    mutate(col_idx = row_number()) %>%
    select(col_idx, survtype) %>%
    filter(survtype == "DHS") %>%
    mutate(row_idx = row_number(),
           x=1)

  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$observations$full_obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)

  mics_join <- mf$observations$full_obs %>%
    mutate(col_idx = row_number()) %>%
    select(col_idx, survtype) %>%
    filter(survtype == "MICS") %>%
    mutate(row_idx = row_number(),
           x=1)

  X_extract_mics <- spMatrix(nrow(mics_join), nrow(mf$observations$full_obs), i=mics_join$row_idx, j=mics_join$col_idx, x=mics_join$x)

  X_extract <- list()
  X_extract$X_extract_ais <- X_extract_ais
  X_extract$X_extract_dhs <- X_extract_dhs
  X_extract$X_extract_mics <- X_extract_mics

  R <- list()
  R$R_spatial <- make_adjacency_matrix(areas, naomi_level)
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips_dhs), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  R$R_country <- as(diag(1, 1), "dgTMatrix")
  R$R_spatial_iid <- as(diag(1, length(unique(mf$mf_model$area_id))), "dgTMatrix")

  mf$mics_toggle <- 0

  if(nrow(mics_join)) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$observations$full_obs %>% filter(survtype == "MICS"))
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% filter(survtype == "MICS"))
    R$R_tips_mics <- make_rw_structure_matrix(ncol(Z$Z_tips_mics), 1, adjust_diagonal = TRUE)
    mf$mics_toggle <- 1
  }

  mf$R <- R
  mf$Z <- Z
  mf$X_extract <- X_extract
  mf$M_naomi_obs <- M_naomi_obs
  mf$M_full_obs <- M_full_obs

  return(mf)
}

extend_populations <- function(population, areas) {

  pop <- population %>%
    dplyr::mutate(period = as.numeric(naomi::year_labels(naomi::calendar_quarter_to_quarter_id(calendar_quarter)))) %>%
    dplyr::select(area_id, period, age_group, sex, population)

  attr(pop, "problems") <- NULL

  pop <- dfertility::area_populations(pop, naomi::spread_areas(areas), 2020) %>%
    dplyr::filter(sex == "female") %>%
    dplyr::ungroup() %>%
    dplyr::select(-sex) %>%
    dplyr::rename(value = population)

  pop <- tidyr::crossing(area_id = unique(pop$area_id),
                                age_group = unique(pop$age_group),
                                period = 1995:2020
  ) %>%
    dplyr::left_join(pop, c("area_id", "age_group", "period")) %>%
    dplyr::group_by(area_id, age_group) %>%
    dplyr::mutate(value = exp(zoo::na.approx(log(value), period, na.rm = FALSE))) %>%
    dplyr::group_by(area_id, age_group) %>%
    tidyr::fill(value, .direction="updown")

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
                              model_level,
                              project = FALSE) {

  areas_long <- areas %>%
    st_drop_geometry() %>%
    mutate(iso3 = iso3_c)

  max_year <- project

  population <- extend_populations(population, areas) %>%
    dplyr::rename(population = value) %>%
    dplyr::mutate(iso3 = iso3_c)

  area_tree <- create_areas(area_merged = areas)
  area_aggregation <- create_area_aggregation(areas$area_id[areas$area_level == model_level], area_tree)

  mf_model <- tidyr::crossing(period = 1995:max_year,
                       age_group = unique(asfr$age_group),
                       area_id = unique(area_aggregation$model_area_id)) %>%
    left_join(
      population %>%
        dplyr::select(iso3, area_id, period, age_group, population),
      by = c("period", "age_group", "area_id")
    ) %>%
    mutate(area_id = factor(area_id),
           age_group = factor(age_group, levels = unique(asfr$age_group)),
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
           age_group = factor(age_group, levels(mf_model$age_group)),
           area_id = factor(area_id, levels(mf_model$area_id)),
           period = factor(period, levels(mf_model$period)),
           spike_2000 = ifelse(period == 2000, 1, 0),
           spike_1999 = ifelse(period == 1999, 1, 0),
           spike_2001 = ifelse(period == 2001, 1, 0)
    )

  mf <- list()
  mf$mf_model <- mf_model
  mf$district$obs <- obs
  mf$mics_toggle <- 0
  mf$out_toggle <- 0

  ## Outputs

  age_aggregation <- data.frame(
    "age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group,
    "model_age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
  ) %>%
    bind_rows(data.frame(
      "age_group" = "15-49",
      "model_age_group" = filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
    ))

  asfr_out <- tidyr::crossing(
    area_id = area_aggregation$area_id,
    age_group = unique(mf_model$age_group),
    period = unique(mf_model$period),
    variable = "asfr"
  ) %>%
    arrange(variable, area_id, age_group, period) %>%
    mutate(out_idx = row_number()) %>%
    droplevels()

  asfr_join_out <- tidyr::crossing(area_aggregation,
                            age_group = unique(mf_model$age_group),
                            period = unique(mf_model$period)) %>%
    full_join(mf_model %>%
                dplyr::select(area_id, age_group, period, idx),
              by = c("model_area_id" = "area_id", "age_group", "period")
    ) %>%
    full_join(asfr_out, by = c("area_id", "age_group", "period")) %>%
    mutate(x=1) %>%
    filter(!is.na(model_area_id))

  tfr_out <- tidyr::crossing(
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

  tfr_join_out <- tidyr::crossing(area_id = area_aggregation$area_id,
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


  if(!is.null(mics_asfr)) {

    # mf_mics <- tidyr::crossing(area_id = filter(areas_long, area_level == 1)$area_id,
    mf_mics <- tidyr::crossing(area_id = as.character(unique(mics_asfr$area_id)),
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
             idx =factor(idx, levels(mf_mics$idx)),
             spike_2000 = ifelse(period == 2000, 1, 0),
             spike_1999 = ifelse(period == 1999, 1, 0),
             spike_2001 = ifelse(period == 2001, 1, 0)
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
  Z$X_urban_dummy <- model.matrix(~0 + urban, mf$mf_model)

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
  R$R_spatial <- make_adjacency_matrix(areas, model_level)
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  R$R_country <- as(diag(1, 1), "dgTMatrix")
  R$R_spatial_iid <- as(diag(1, length(unique(mf$mf_model$area_id))), "dgTMatrix")

  if(!is.null(mics_asfr)) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs)
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    # X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
    R$R_tips_mics <- make_rw_structure_matrix(ncol(Z$Z_tips_mics), 1, adjust_diagonal = TRUE)

    mf$mics$M_obs_mics <- M_obs_mics
  }

  mf$R <- R
  mf$Z <- Z
  mf$X_extract <- X_extract

  mf$district$M_obs <- M_obs


  return(mf)
}

tmb_outputs <- function(fit, mf, areas) {

  areas_long <- areas %>% sf::st_drop_geometry()

  asfr_qtls <- apply(fit$sample$lambda_out, 1, quantile, c(0.025, 0.5, 0.975))
  tfr_qtls <- apply(fit$sample$tfr_out, 1, quantile, c(0.025, 0.5, 0.975))

  tmb_results <- mf$out$mf_out %>%
    dplyr::mutate(lower = c(asfr_qtls[1,], tfr_qtls[1,]),
                  median = c(asfr_qtls[2,], tfr_qtls[2,]),
                  upper = c(asfr_qtls[3,], tfr_qtls[3,]),
                  source = "tmb") %>%
    utils::type.convert() %>%
    dplyr::left_join(areas_long)

  return(tmb_results)

}
