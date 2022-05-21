#' Aggregate district level populations
#' @description Aggregate district level populations to all levels of area hierarchy
#' @param population Age/sex/space/time stratified population
#' @param areas_wide Area hierarchy in wide format
#' @export

area_populations <- function(population, areas_wide, project, naomi_level) {

  population <- population %>%
    dplyr::filter(sex == "female", period <= project)

  col_nam <- c(paste0("area_id", 0:naomi_level),
    paste0("area_name", 0:naomi_level)
  )

  areas_wide <- areas_wide %>%
    dplyr::select(all_of(col_nam)) %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct() %>%
    dplyr::mutate(area_id = .data[[paste0("area_id", naomi_level)]])

  base_area_pop <- dplyr::left_join(
    areas_wide, population, by = c("area_id"))

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

  R_mat <- methods::as(R_mat, "dgCMatrix")

  return(R_mat)

}

extract_model_level <- function(df_list, model_level) {
  df_list %>%
    dplyr::filter(area_level == model_level) %>%
    tidyr::separate(area_id, into=c("iso3", NA), sep=3, remove=FALSE)
}

make_adjacency_matrix <- function(areas, model_level) {

  if(length(model_level) > 1) {
   filtered_areas <- Map(extract_model_level, areas, model_level) %>%
     dplyr::bind_rows()
  } else {
    filtered_areas <- areas %>%
      dplyr::bind_rows() %>%
      dplyr::filter(area_level == model_level)
  }

  sh <- filtered_areas %>%
    dplyr::group_by(iso3) %>%
    dplyr::arrange(area_sort_order, .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(area_idx = row_number())

  #' Neighbor list
  nb <- sh %>%
    # st_as_sf() %>%
    # methods::as("Spatial") %>%
    spdep::poly2nb() %>%
    `names<-`(sh$area_idx)

  names(nb) <- sh$area_id


  nb <- lapply(nb, as.integer)
  class(nb) <- "nb"

  adj <- spdep::nb2mat(nb, zero.policy=TRUE, style="B")
  R_spatial <- INLA::inla.scale.model(
    diag(rowSums(adj)) - 0.99*adj,
    constr = list(A = matrix(1, 1, nrow(adj)), e = 0)
  )

  return(R_spatial)
}

#' Make model frames
#' @description Create model frames and aggregation matrices for TMB model.
#' @param iso3 iso3 code for country
#' @param population Age/sex/space stratified population
#' @param asfr ASFRs by district and time
#' @param areas Area hierarchy
#' @param naomi_level Area level to produce estimates at
#' @param project Model will by default produce estimates up to year of last survey. Integer projection year if desired, else FALSE
#' @export
#'

make_model_frames_dev <- function(iso3,
                              population,
                              asfr,
                              areas,
                              naomi_level,
                              project = 2020) {

  areas_long <- areas %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(iso3 = iso3)

  # population <- population %>%
  #   mutate(period = as.numeric(year_labels(calendar_quarter_to_quarter_id(calendar_quarter)))) %>%
  #   select(-calendar_quarter)
  #
  # population <- area_populations(population, spread_areas(areas), project, naomi_level) %>%
  #   dplyr::filter(sex == "female") %>%
  #   ungroup %>%
  #   dplyr::select(-sex)
  #
  # population <- tidyr::crossing(area_id = unique(population$area_id),
  #                        age_group = unique(population$age_group),
  #                        period = 1995:2020
  # ) %>%
  #   left_join(population, c("area_id", "age_group", "period")) %>%
  #   group_by(area_id, age_group) %>%
  #   mutate(population = exp(zoo::na.approx(log(population), period, na.rm = FALSE))) %>%
  #   tidyr::fill(population, .direction="updown") %>%
  #   left_join(
  #     areas_long %>% dplyr::select(area_id, iso3),
  #     by = "area_id"
  #   ) %>%
  #   mutate(population = ifelse(is.nan(population), 0, population))

  areas <- dplyr::filter(areas, area_level <= naomi_level)

  area_tree <- naomi::create_areas(area_merged = areas)
  area_aggregation <- naomi::create_area_aggregation(areas$area_id[areas$area_level == naomi_level], area_tree)

  mf_model <- tidyr::crossing(period = 1995:project,
                       age_group = unique(asfr$age_group),
                       area_id = unique(area_aggregation$model_area_id)) %>%
    dplyr::left_join(
      population %>%
        dplyr::select(area_id, period, age_group, population),
      by = c("period", "age_group", "area_id")
    ) %>%
    dplyr::mutate(area_id = factor(area_id),
           age_group = factor(age_group, levels = unique(asfr$age_group)),
           period = factor(period),
           urban = ifelse(area_id %in% c(
             dplyr::filter(areas_long, parent_area_id == "ETH_1_10")$area_id,
             dplyr::filter(areas_long, str_detect(area_name, "Town"))$area_id,
             dplyr::filter(areas_long, area_name %in% c("Harari", "Fafen (Jijiga)", "Erer"))$area_id),
             1, 0)
    ) %>%
    dplyr::arrange(period, area_id, age_group) %>%
    dplyr::mutate(idx = factor(row_number()),
           id.period = group_indices(., period)-1,
           id.interaction1 = factor(group_indices(., age_group, period)),
           id.interaction2 = factor(group_indices(., period, area_id)),
           id.interaction3 = factor(group_indices(., age_group, area_id))
           # id.omega1 = factor(group_indices(., age_group, iso3)),
           # id.omega2 = factor(group_indices(., period, iso3))
    )

  obs <- asfr %>%
    dplyr::left_join(areas_long) %>%
    dplyr::mutate(period = factor(period, levels(mf_model$period))) %>%
    dplyr::filter(!is.na(births)) %>%
    dplyr::select(survey_id, survtype, area_id, area_level, period, age_group, tips, births, pys) %>%
    # left_join(full_frame) %>%
    dplyr::mutate(tips_dummy = as.integer(tips > 5),
           tips_dummy_10 = as.integer(tips == 10),
           tips_dummy_9_11 = as.integer(tips %in% c(9, 11)),
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
    dplyr::filter(area_level == naomi_level) %>%
    dplyr::mutate(
           area_id = factor(area_id, levels(mf_model$area_id))
    ) %>%
    dplyr::left_join(mf_model %>% dplyr::select(area_id, age_group, period, idx))

  aggregate_mf <- tidyr::crossing(area_id = area_aggregation$area_id,
                         period = unique(mf_model$period),
                         age_group = unique(mf_model$age_group)
  ) %>%
    dplyr::arrange(area_id, age_group, period) %>%
    dplyr::mutate(idx_row = factor(row_number()))

  join <- aggregate_mf %>%
    dplyr::left_join(area_aggregation, by = "area_id") %>%
    dplyr::left_join(mf_model, by = c("age_group", "period", "model_area_id" = "area_id")) %>%
    dplyr::mutate(x = 1) %>%
    utils::type.convert()

  full_obs <- obs %>%
    dplyr::left_join(aggregate_mf, by=c("area_id", "age_group", "period")) %>%
    dplyr::rename(idx = idx_row)

  A_full_obs <- Matrix::sparseMatrix(
    i = join$idx_row,
    j = join$idx,
    x = join$x,
    dims = c(max(join$idx_row), nrow(mf_model)),
    use.last.ij = TRUE
  )

  mf <- list()
  mf$mf_model <- mf_model
  mf$observations$naomi_level_obs <- naomi_level_obs
  mf$observations$full_obs <- full_obs
  mf$observations$A_full_obs <- A_full_obs

  ## Outputs

  age_aggregation <- data.frame(
    "age_group" = dplyr::filter(naomi::get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group,
    "model_age_group" = dplyr::filter(naomi::get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
  ) %>%
    dplyr::bind_rows(data.frame(
      "age_group" = "Y015_Y049",
      "model_age_group" = dplyr::filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
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
  #   dplyr::filter(!is.na(model_area_id))

  tfr_out <- tidyr::crossing(
    area_id = area_aggregation$area_id,
    age_group = "Y015_Y049",
    period = unique(mf_model$period),
    variable = "tfr"
  ) %>%
    dplyr::arrange(variable, area_id, age_group, period) %>%
    dplyr::mutate(idx_out = dplyr::row_number()) %>%
    droplevels()

  tfr_join_out <- tidyr::crossing(area_id = area_aggregation$area_id,
                           age_aggregation %>% dplyr::filter(age_group == "Y015_Y049"),
                           period = unique(mf_model$period)) %>%
    dplyr::full_join(
      aggregate_mf %>%
        dplyr::select(area_id, age_group, period, idx_row) %>%
        unique,
      by = c("area_id", "model_age_group" = "age_group", "period")
    ) %>%
    dplyr::rename(idx_col = idx_row) %>%
    dplyr::arrange(period, area_id, age_group) %>%
    dplyr::full_join(tfr_out, by = c("area_id", "age_group", "period")) %>%
    dplyr::mutate(x = 5)

  # A_asfr_out <- spMatrix(nrow(aggregate_mf), nrow(mf_model), asfr_join_out$out_idx, as.integer(asfr_join_out$idx), asfr_join_out$x)
  A_tfr_out <- Matrix::spMatrix(nrow(tfr_out), nrow(aggregate_mf), tfr_join_out$idx_out, as.integer(tfr_join_out$idx_col), tfr_join_out$x)

  mf_out <- aggregate_mf %>%
    dplyr::mutate(variable = "asfr") %>%
    dplyr::bind_rows(tfr_out) %>%
    dplyr::select(-idx_out)

  mf$out$mf_out <- mf_out
  # mf$out$A_asfr_out <- A_asfr_out
  mf$out$A_tfr_out <- A_tfr_out

  M_naomi_obs <- Matrix::sparse.model.matrix(~0 + idx, mf$observations$naomi_level_obs)
  M_full_obs <- Matrix::sparse.model.matrix(~0 + idx, mf$observations$full_obs)

  Z <- list()
  Z$Z_spatial <- Matrix::sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- Matrix::sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- Matrix::sparse.model.matrix(~0 + period, mf$mf_model)
  Z$Z_country <- methods::as(matrix(rep(1, nrow(mf$mf_model)), ncol = 1), "dgCMatrix")

  # Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$observations$full_obs)
  Z$Z_tips_dhs <- Matrix::sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% dplyr::filter(survtype == "DHS"))
  Z$Z_tips_ais <- Matrix::sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% dplyr::filter(survtype %in% c("AIS", "MIS")))
  Z$X_tips_dummy <- Matrix::model.matrix(~0 + tips_dummy, mf$observations$full_obs %>% dplyr::filter(survtype == "DHS"))
  Z$X_tips_dummy_10 <- Matrix::model.matrix(~0 + tips_dummy_10, mf$observations$full_obs %>% dplyr::filter(survtype == "DHS"))
  Z$X_tips_dummy_9_11 <- Matrix::model.matrix(~0 + tips_dummy_9_11, mf$observations$full_obs %>% dplyr::filter(survtype == "DHS"))
  Z$X_urban_dummy <- Matrix::model.matrix(~0 + urban, mf$mf_model)
  Z$X_period <- methods::as(Matrix::as.matrix(mf_model$id.period), "dgTMatrix")

  ais_join <- mf$observations$full_obs %>%
    dplyr::mutate(col_idx = dplyr::row_number()) %>%
    dplyr::select(col_idx, survtype) %>%
    dplyr::filter(survtype %in% c("AIS", "MIS")) %>%
    dplyr::mutate(row_idx = dplyr::row_number(),
           x=1)

  X_extract_ais <- Matrix::spMatrix(
    nrow(ais_join),
    nrow(mf$observations$full_obs),
    i = ais_join$row_idx,
    j = ais_join$col_idx,
    x = ais_join$x
  )

  dhs_join <- mf$observations$full_obs %>%
    dplyr::mutate(col_idx = dplyr::row_number()) %>%
    dplyr::select(col_idx, survtype) %>%
    dplyr::filter(survtype == "DHS") %>%
    dplyr::mutate(row_idx = dplyr::row_number(),
           x=1)

  X_extract_dhs <- Matrix::spMatrix(
    nrow(dhs_join),
    nrow(mf$observations$full_obs),
    i = dhs_join$row_idx,
    j = dhs_join$col_idx, x=dhs_join$x
  )

  phia_join <- mf$observations$full_obs %>%
    dplyr::mutate(col_idx = dplyr::row_number()) %>%
    dplyr::select(col_idx, survtype) %>%
    dplyr::filter(survtype == "PHIA") %>%
    dplyr::mutate(row_idx = dplyr::row_number(),
           x=1)

  X_extract_phia <- Matrix::spMatrix(
    nrow(phia_join),
    nrow(mf$observations$full_obs),
    i = phia_join$row_idx,
    j = phia_join$col_idx, x=phia_join$x
  )

  mics_join <- mf$observations$full_obs %>%
    dplyr::mutate(col_idx = dplyr::row_number()) %>%
    dplyr::select(col_idx, survtype) %>%
    dplyr::filter(survtype == "MICS") %>%
    dplyr::mutate(row_idx = dplyr::row_number(),
           x=1)

  X_extract_mics <- Matrix::spMatrix(
    nrow(mics_join),
    nrow(mf$observations$full_obs),
    i = mics_join$row_idx,
    j = mics_join$col_idx,
    x = mics_join$x
  )

  X_extract <- list()
  X_extract$X_extract_ais <- X_extract_ais
  X_extract$X_extract_dhs <- X_extract_dhs
  X_extract$X_extract_mics <- X_extract_mics
  X_extract$X_extract_phia <- X_extract_phia

  R <- list()
  R$R_spatial <- make_adjacency_matrix(areas, naomi_level)
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips_dhs), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  R$R_country <- methods::as(diag(1, 1), "dgTMatrix")
  R$R_spatial_iid <- methods::as(diag(1, length(unique(mf$mf_model$area_id))), "dgTMatrix")

  mf$mics_toggle <- 0

  if(nrow(mics_join)) {
    M_obs_mics <- Matrix::sparse.model.matrix(~0 + idx, mf$observations$full_obs %>% dplyr::filter(survtype == "MICS"))
    Z$Z_tips_mics <- Matrix::sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% dplyr::filter(survtype == "MICS"))
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

#' Make model frames for batch running
#' @description Create model frames and aggregation matrices for TMB model.
#' @param population Age/sex/space stratified population
#' @param asfr ASFRs by district and time
#' @param areas_list List of area files
#' @param lvl Dataframe of district and province levels
#' @param project Model will by default produce estimates up to year of last survey. Integer projection year if desired, else FALSE
#' @export

make_model_frames_batch <- function(lvl_map,
                                    population,
                                    asfr,
                                    areas_list,
                                    project = 2020) {

  fertility_fit_level <- lvl_map$fertility_fit_level

  areas <- areas_list %>%
    dplyr::bind_rows()
  #
  # population <- population %>%
  #   bind_rows()

  area_aggregation <- Map(function(areas, fertility_fit_level) {
    areas <- dplyr::filter(areas, area_level <= fertility_fit_level)
    area_tree <- naomi::create_areas(area_merged = areas)
    area_aggregation <- naomi::create_area_aggregation(areas$area_id[areas$area_level == fertility_fit_level], area_tree) %>%
      dplyr::left_join(
        areas %>%
          dplyr::select(area_id, area_sort_order) %>%
          sf::st_drop_geometry(),
        by = c("model_area_id" = "area_id")
      )
  }, areas_list, fertility_fit_level) %>%
    dplyr::bind_rows(.id = "iso3") %>%
    # left_join(areas %>% select(iso3, area_id, area_sort_order) %>% st_drop_geometry(), by=c("model_area_id" = "area_id")) %>%
   dplyr::group_by(iso3) %>%
   dplyr::arrange(area_sort_order, .by_group = TRUE) %>%
   dplyr::ungroup()

  # asfr <- asfr %>%
  #   bind_rows()

  # spatial_interaction_merge <- areas_list %>%
  #   lapply(spread_areas) %>%
  #   Map(function(areas, interaction_level) {
  #     df <- areas %>%
  #       select(area_id, starts_with("area_id") & ends_with(as.character(interaction_level))) %>%
  #       st_drop_geometry()
  #
  #     colnames(df) <- c("area_id", "spatial_area_id")
  #
  #     df
  #   }, ., interaction_level) %>%
  #   bind_rows()


  mf_model <- tidyr::crossing(period = 1995:project,
                       age_group = unique(asfr$age_group),
                       area_id = unique(area_aggregation$model_area_id)) %>%
    dplyr::left_join(
      population %>%
        dplyr::select(iso3, area_id, period, age_group, population),
      by = c("period", "age_group", "area_id")
    ) %>%
    # left_join(spatial_interaction_merge) %>%
    dplyr::mutate(area_id = factor(area_id, levels = unique(area_aggregation$model_area_id)),
           age_group = factor(age_group, levels = unique(asfr$age_group)),
           period = factor(period)
           # urban = ifelse(area_id %in% c(
           #   dplyr::filter(areas_long, parent_area_id == "ETH_1_10")$area_id,
           #   dplyr::filter(areas_long, str_detect(area_name, "Town"))$area_id,
           #   dplyr::filter(areas_long, area_name %in% c("Harari", "Fafen (Jijiga)", "Erer"))$area_id),
           #   1, 0)
    ) %>%
    dplyr::arrange(period, area_id, age_group) %>%
    dplyr::mutate(idx = factor(dplyr::row_number()),
           id.period = dplyr::group_indices(., period)-1,
           id.interaction1 = factor(dplyr::group_indices(., age_group, period, iso3)),
           id.interaction2 = factor(dplyr::group_indices(., period, area_id)),
           id.interaction3 = factor(dplyr::group_indices(., age_group, area_id)),
           id.omega1 = factor(dplyr::group_indices(., age_group, iso3)),
           id.omega2 = factor(dplyr::group_indices(., period, iso3))
    )

  # if (any(is.na(mf_model$population))) {
  #   missing_population <- mf_model %>%
  #     dplyr::filter(is.na(population)) %>%
  #     select(period, age_group, area_id)
  #   stop("Age, time, space combinations without populations:\n",
  #          paste0(utils::capture.output(missing_population), collapse = "\n"))
  # }

  obs <- asfr %>%
    dplyr::left_join(areas %>% st_drop_geometry()) %>%
    dplyr::mutate(period = factor(period, levels(mf_model$period))) %>%
    dplyr::filter(!is.na(births)) %>%
    dplyr::select(iso3, survey_id, survtype, area_id, area_level, period, age_group, tips, births, pys) %>%
    # dplyr::left_join(full_frame) %>%
    dplyr::mutate(tips_dummy = as.integer(tips > 5),
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
    dplyr::group_split(iso3) %>%
    Map(function(obs, fertility_fit_level) {
      obs %>%
        dplyr::filter(area_level == fertility_fit_level)
    }, ., fertility_fit_level) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(area_id = factor(area_id, levels(mf_model$area_id))) %>%
    dplyr::left_join(mf_model %>% select(area_id, age_group, period, idx))

  aggregate_mf <- tidyr::crossing(area_id = area_aggregation$area_id,
                         period = unique(mf_model$period),
                         age_group = unique(mf_model$age_group)
  ) %>%
    dplyr::arrange(area_id, age_group, period) %>%
    dplyr::mutate(idx_row = factor(row_number()))

  join <- aggregate_mf %>%
    dplyr::left_join(area_aggregation, by = "area_id") %>%
    dplyr::left_join(mf_model, by=c("iso3", "age_group", "period", "model_area_id" = "area_id")) %>%
    dplyr::mutate(x = 1) %>%
    type.convert()

  full_obs <- obs %>%
    dplyr::left_join(aggregate_mf, by=c("area_id", "age_group", "period")) %>%
    dplyr::rename(idx = idx_row)
  #
  A_full_obs <- Matrix::sparseMatrix(
    i = join$idx_row,
    j = join$idx,
    x = join$x,
    dims = c(max(join$idx_row), nrow(mf_model)),
    use.last.ij = TRUE
  )

  mf <- list()
  mf$mf_model <- mf_model
  mf$observations$naomi_level_obs <- naomi_level_obs
  mf$observations$full_obs <- full_obs
  mf$observations$A_full_obs <- A_full_obs

  ## Outputs

  age_aggregation <- data.frame(
    "age_group" = dplyr::filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group,
    "model_age_group" = dplyr::filter(get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
  ) %>%
    dplyr::bind_rows(data.frame(
      "age_group" = "Y015_Y049",
      "model_age_group" = dplyr::filter(naomi::get_age_groups(), age_group_start %in% 15:45, age_group_span == 5)$age_group
    ))

  # asfr_out <- tidyr::crossing(
  #   area_id = area_aggregation$area_id,
  #   age_group = unique(mf_model$age_group),
  #   period = unique(mf_model$period),
  #   variable = "asfr"
  # ) %>%
  #   dplyr::arrange(variable, area_id, age_group, period) %>%
  #   dplyr::mutate(out_idx = row_number()) %>%
  #   dplyr::droplevels()
  #
  # asfr_join_out <- tidyr::crossing(area_aggregation,
  #                           age_group = unique(mf_model$age_group),
  #                           period = unique(mf_model$period)) %>%
  #   dplyr::full_join(mf_model %>%
  #               dplyr::select(area_id, age_group, period, idx),
  #             by = c("model_area_id" = "area_id", "age_group", "period")
  #   ) %>%
  #   dplyr::full_join(asfr_out, by = c("area_id", "age_group", "period")) %>%
  #   dplyr::mutate(x = 1) %>%
  #   dplyr::filter(!is.na(model_area_id))

  tfr_out <- crossing(
    area_id = area_aggregation$area_id,
    age_group = "Y015_Y049",
    period = unique(mf_model$period),
    variable = "tfr"
  ) %>%
    dplyr::arrange(variable, area_id, age_group, period) %>%
    dplyr::mutate(idx_out = dplyr::row_number()) %>%
    droplevels()

  tfr_join_out <- tidyr::crossing(area_id = area_aggregation$area_id,
                           age_aggregation %>% dplyr::filter(age_group == "Y015_Y049"),
                           period = unique(mf_model$period)) %>%
    dplyr::full_join(
      aggregate_mf %>%
        dplyr::select(area_id, age_group, period, idx_row) %>%
        unique,
      by = c("area_id", "model_age_group" = "age_group", "period")
    ) %>%
    dlpyr::rename(idx_col = idx_row) %>%
    dlpyr::arrange(period, area_id, age_group) %>%
    dlpyr::full_join(tfr_out, by = c("area_id", "age_group", "period")) %>%
    dlpyr::mutate(x = 5)

  # A_asfr_out <- spMatrix(nrow(aggregate_mf), nrow(mf_model), asfr_join_out$out_idx, as.integer(asfr_join_out$idx), asfr_join_out$x)
  A_tfr_out <- Matrix::spMatrix(
    nrow(tfr_out),
    nrow(aggregate_mf),
    tfr_join_out$idx_out,
    as.integer(tfr_join_out$idx_col),
    tfr_join_out$x
  )

  mf_out <- aggregate_mf %>%
    dplyr::mutate(variable = "asfr") %>%
    dplyr::bind_rows(tfr_out) %>%
    dplyr::select(-idx_out)

  mf$out$mf_out <- mf_out
  # mf$out$A_asfr_out <- A_asfr_out
  mf$out$A_tfr_out <- A_tfr_out

  M_naomi_level_obs <- sparse.model.matrix(~0 + idx, mf$observations$naomi_level_obs)
  M_full_obs <- sparse.model.matrix(~0 + idx, mf$observations$full_obs)

  Z <- list()
  Z$Z_spatial <- Matrix::sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- Matrix::sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- Matrix::sparse.model.matrix(~0 + period, mf$mf_model)
  Z$Z_country <- Matrix::sparse.model.matrix(~0 + iso3, mf$mf_model)

  # Z$Z_tips <- Matrix::sparse.model.matrix(~0 + tips_f, mf$observations$full_obs)
  Z$Z_tips_dhs <- Matrix::sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% dplyr::filter(survtype == "DHS"))
  Z$Z_tips_ais <- Matrix::sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% dplyr::filter(survtype %in% c("AIS", "MIS")))
  Z$X_tips_dummy <- stats::model.matrix(~0 + tips_dummy, mf$observations$full_obs %>% dplyr::filter(survtype == "DHS"))
  # Z$X_urban_dummy <- model.matrix(~0 + urban, mf$mf_model)
  Z$X_period <- methods::as(Matrix::as.matrix(mf_model$id.period), "dgTMatrix")

  ais_join <- mf$observations$full_obs %>%
    dplyr::mutate(col_idx = dplyr::row_number()) %>%
    dplyr::select(col_idx, survtype) %>%
    dplyr::filter(survtype %in% c("AIS", "MIS")) %>%
    dplyr::mutate(row_idx = row_number(),
           x=1)

  X_extract_ais <-  Matrix::spMatrix(
    nrow(ais_join),
    nrow(mf$observations$full_obs),
    i = ais_join$row_idx,
    j = ais_join$col_idx,
    x = ais_join$x
  )

  dhs_join <- mf$observations$full_obs %>%
   dplyr::mutate(col_idx = dplyr::row_number()) %>%
   dplyr::select(col_idx, survtype) %>%
   dplyr::filter(survtype == "DHS") %>%
   dplyr::mutate(row_idx = dplyr::row_number(), x = 1)

  X_extract_dhs <- Matrix::spMatrix(
    nrow(dhs_join),
    nrow(mf$observations$full_obs),
    i = dhs_join$row_idx,
    j = dhs_join$col_idx,
    x = dhs_join$x
  )

  mics_join <- mf$observations$full_obs %>%
    dplyr::mutate(col_idx = dplyr::row_number()) %>%
    dplyr::select(col_idx, survtype) %>%
    dplyr::filter(survtype == "MICS") %>%
    dplyr::mutate(row_idx = dplyr::row_number(), x = 1)

  X_extract_mics <- Matrix::spMatrix(
    nrow(mics_join),
    nrow(mf$observations$full_obs),
    i = mics_join$row_idx,
    j = mics_join$col_idx,
    x = mics_join$x
  )

  X_extract <- list()
  X_extract$X_extract_ais <- X_extract_ais
  X_extract$X_extract_dhs <- X_extract_dhs
  X_extract$X_extract_mics <- X_extract_mics

  R <- list()
  R$R_spatial <-  make_adjacency_matrix(areas_list, fertility_fit_level)
  # R$R_spatial_interaction <- make_adjacency_matrix(areas_list, interaction_level)
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips_dhs), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 1, adjust_diagonal = TRUE)
  R$R_country <- methods::as(diag(1, length(areas_list)), "dgTMatrix")
  R$R_spatial_iid <- methods::as(diag(1, length(unique(mf$mf_model$area_id))), "dgTMatrix")

  mf$mics_toggle <- 0

  if(nrow(mics_join)) {
    M_obs_mics <- Matrix::sparse.model.matrix(~0 + idx, mf$observations$full_obs %>% dplyr::filter(survtype == "MICS"))
    Z$Z_tips_mics <- Matrix::sparse.model.matrix(~0 + tips_f, mf$observations$full_obs %>% dplyr::filter(survtype == "MICS"))
    R$R_tips_mics <- make_rw_structure_matrix(ncol(Z$Z_tips_mics), 1, adjust_diagonal = TRUE)
    mf$mics_toggle <- 1
  }

  mf$R <- R
  mf$Z <- Z
  mf$X_extract <- X_extract
  mf$M_naomi_level_obs <- M_naomi_level_obs
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
    tidyr::fill(value, .direction = "updown")

}

validate_model_frame <- function(mf, areas) {

  population_check <- mf$mf_model %>%
    dplyr::filter(is.na(population))

  if(nrow(population_check))
    stop("Areas in model frame have no population")

  observation_area_id <- unique(as.character(mf$observations$full_obs$area_id))
  area_check <- observation_area_id[!observation_area_id %in% areas$area_id]

  if(!rlang::is_empty(area_check))
    stop(paste0(
      "Area IDs in observation dataset are not found in area file\n\n",
      area_check
      )
    )

}


tmb_outputs <- function(fit, mf, areas) {

  areas_long <- sf::st_drop_geometry(areas)

  asfr_qtls <- apply(fit$sample$lambda_out, 1, stats::quantile, c(0.025, 0.5, 0.975))
  tfr_qtls <- apply(fit$sample$tfr_out, 1, stats::quantile, c(0.025, 0.5, 0.975))

  tmb_results <- mf$out$mf_out %>%
    dplyr::mutate(lower = c(asfr_qtls[1,], tfr_qtls[1,]),
                  median = c(asfr_qtls[2,], tfr_qtls[2,]),
                  upper = c(asfr_qtls[3,], tfr_qtls[3,]),
                  source = "tmb") %>%
    utils::type.convert() %>%
    dplyr::left_join(areas_long)

  return(tmb_results)

}
