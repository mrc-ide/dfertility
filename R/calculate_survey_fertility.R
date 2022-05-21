get_fertility_surveys <- function(surveys) {

  ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)

  ird <- ird %>%
    dplyr::mutate(path = unlist(get_datasets(., clear_cache = TRUE))) %>%
    dplyr::bind_rows()

  ir <- lapply(ird$path, readRDS) %>%
    lapply(function(x) {class(x) <- "data.frame"
    return(x)}) %>%
    Map(function(ir, surveys) {
      dplyr::mutate(ir,
             surveyid = surveys$survey_id,
             country = surveys$CountryName,
             survyear = surveys$SurveyYear,
             survtype = surveys$SurveyType)
    }, ., dplyr::group_split(surveys, SurveyId))

  out <- list()
  out$ir <- ir
  out$nrow_ir <- lapply(ir, nrow)
  out


}

assign_tips <- function(ir, single_tips) {

  survey_type <- ir %>%
    lapply("[", "survtype") %>%
    lapply(unique) %>%
    dplyr::bind_rows()

  if(!single_tips)
    # tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[survey_type$survtype]
    tips_surv <- lapply(survey_type$survtype, function(x) {
      if(x=="DHS")
        c(0,15)
      else
        c(0,5)
    })
  else
    # tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[survey_type$survtype]
    tips_surv <- lapply(survey_type$survtype, function(x) {
      if(x=="DHS")
        c(0:15)
      else
        c(0:5)
    })

  tips_surv
}

map_ir_to_areas <- function(ir, cluster_list, single_tips = TRUE) {

  mc.cores <- if(.Platform$OS.type == "windows") 1 else parallel::detectCores()

  ir <- Map(ir_by_area, ir, cluster_list[names(ir)], n = 1:length(ir), total = length(ir)) %>%
    unlist(recursive = FALSE)

  tips_surv <- assign_tips(ir, single_tips)

  dat <- list()
  dat$ir <- ir
  dat$tips_surv <- tips_surv

  dat
}

#' Aggregate district level cluster dataframe
#' @description Aggregate a dataframe of survey clusters geolocated to an administrative to any higher administrative level
#' @param clusters Dataframe of geolocated survey clusters
#' @param areas_wide Area hierarchy in wide format
#' @param area_level Integer referring to administrative level, where 0 is national, 1 is provincial etc.
#' @return Dataframe list of clusters by administrative area, one list item per survey
#' @seealso \code{\link{get_areas}}
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export

assign_cluster_area <- function(clusters, areas_wide, area_level) {

  areas <- clusters %>%
    dplyr::rename(area_id = geoloc_area_id) %>%
    dplyr::left_join(areas_wide %>% dplyr::select(area_id, paste0("area_id", area_level))) %>%
    dplyr::select(-area_id) %>%
    dplyr::rename(area_id = paste0("area_id", area_level))

  area_list <- areas %>%
    dplyr::group_by(survey_id) %>%
    dplyr::group_split(keep = TRUE)

  names(area_list) <- area_list %>%
    lapply("[", "survey_id") %>%
    lapply(unique) %>%
    dplyr::bind_rows %>%
    .$survey_id

  return(area_list)

}


#' Join individual recode survey datasets by area
#' @param ir Individual recode survey dataset
#' @param area_list List of areas
ir_by_area <- function(ir, area_list, n, total) {

  print(paste(n, "of", total))

  ir_int <- ir %>%
    left_join(area_list, by=c("v001" = "cluster_id")) %>%
    mutate(survey_id = factor(survey_id),
           survtype = factor(survtype),
           survyear = factor(survyear),
           area_id = factor(area_id)) %>%
    dplyr::filter(!is.na(area_id)) %>%
    group_split(area_id)

  return(ir_int)

}

# #' Join individual recode survey datasets by cluster to area IDs
# #' @param surveys Survey dataframe as created by \code{rdhs::dhs_surveys()}
# #' @param cluster_areas Dataframe of clusters assigned to areas. Output of \code{\link{cluster_areas}}.
# #' @param single_tips desc this
# #' @return Nested list: \code{df$ir} and \code{df$tips_surv}. \code{df$ir} has one list item per survey/area combination and \code{df$tips_surv} has the corroborating TIPS value.
# #' @export
# #' @seealso \code{\link{cluster_areas}}
#
# clusters_to_surveys <- function(iso3, surveys, cluster_areas, level, single_tips = TRUE) {
#
#   ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
#
#   ird <- ird %>%
#     mutate(path = unlist(get_datasets(.))) %>%
#     bind_rows()
#
#   ir <- lapply(ird$path, readRDS) %>%
#     lapply(function(x) {class(x) <- "data.frame"
#     return(x)}) %>%
#     Map(function(ir, surveys) {
#       mutate(ir,
#              surveyid = surveys$survey_id,
#              country = surveys$CountryName,
#              survyear = surveys$SurveyYear,
#              survtype = surveys$SurveyType)
#     }, ., group_split(surveys, SurveyId))
#
#
#   ## I think this mess was necessary because the order of cluster_areas was not always the same as the order for ir.. Double check. Otherwise the simple names(ir) <- names(cluster_areas) is sufficient
#   # names(ir) <- ir %>%
#   #   lapply("[", "surveyid") %>%
#   #   lapply(unique) %>%
#   #   bind_rows %>%
#   #   left_join(clusters %>%
#   #               select(survey_id, DHS_survey_id) %>%
#   #               unique,
#   #             by=c("surveyid" = "DHS_survey_id")) %>%
#   #   select(-surveyid) %>%
#   #   .$survey_id
#
#   names(ir) <- names(cluster_areas)
#
#   if(iso3 == "ETH") {
#
#     cols_edit <- c("v008", "v011", "b3_01", "b3_02", "b3_03", "b3_04", "b3_05", "b3_06", "b3_07", "b3_08", "b3_09", "b3_10", "b3_11", "b3_12", "b3_13", "b3_14", "b3_15", "b3_16", "b3_17", "b3_18", "b3_19", "b3_20")
#
#     ir <- ir %>%
#       lapply(function(x) x %>% mutate_at(.vars = cols_edit, .funs = ~(.+92)))
#
#   }
#
#
#   if(level > 0) {
#
#     ir <- Map(ir_by_area, ir, cluster_areas[names(ir)], n=1:length(ir), total=length(ir)) %>%
#       unlist(recursive = FALSE)
#
#     survey_type <- ir %>%
#       lapply("[", "survtype") %>%
#       lapply(unique) %>%
#       bind_rows
#
#     if(!single_tips)
#       tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[survey_type$survtype]
#     else
#       tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[survey_type$survtype]
#
#   } else {
#
#     ir <- lapply(ir, function(x) {
#       mutate(x, area_id = unique(surveys$iso3))
#     })
#
#     if(!single_tips)
#       tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[surveys$SurveyType]
#     else
#       tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[surveys$SurveyType]
#
#   }
#
#   dat <- list()
#   dat$ir <- ir
#   dat$tips_surv <- tips_surv
#
#   return(dat)
#
# }

#' Calculate fertility rates from MICS data
#' @export
calculate_mics_fertility <- function(iso3, mics_wm, mics_births_to_women) {

  mc.cores <- if(.Platform$OS.type == "windows") 1 else parallel::detectCores()

  mics_wm_asfr <- mics_wm %>%
    type.convert() %>%
    mutate(survey_id = factor(survey_id),
           area_id = factor(area_id)) %>%
    arrange(survey_id) %>%
    group_by(survey_id) %>%
    group_split()

  mics_births_asfr <- mics_births_to_women %>%
    left_join(mics_wm) %>%
    mutate(survey_id = factor(survey_id),
           area_id = factor(area_id)) %>%
    arrange(survey_id) %>%
    group_by(survey_id) %>%
    group_split()

  message("Calculating district-level MICS ASFR")

  # For model:
  mics_asfr <- Map(calc_asfr, mics_wm_asfr,
                   by = list(~area_id + survey_id),
                   tips = list(c(0:15)),
                   agegr= list(3:10*5),
                   period = list(1995:2019),
                   clusters = list(~cluster),
                   strata = list(NULL),
                   id = list("unique_id"),
                   dob = list("wdob"),
                   intv = list("doi"),
                   weight = list("weight"),
                   varmethod = list("none"),
                   bhdata = mics_births_asfr,
                   bvars = list("cdob"),
                   counts = TRUE) %>%
    bind_rows %>%
    type.convert() %>%
    separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
    dplyr::filter(period <= survyear) %>%
    # rename(age_group = agegr) %>%
    mutate(survtype = "MICS",
           iso3 = iso3
    ) %>%
    left_join(get_age_groups() %>% select(age_group, age_group_label), by=c("agegr" = "age_group_label")) %>%
    select(-agegr)

  message("Calculating aggregate MICS fertility rates")
  # For plotting:
  mics_asfr_plot <- Map(calc_asfr, mics_wm_asfr,
                        by = list(~area_id + survey_id),
                        tips = list(c(0,15)),
                        agegr= list(3:10*5),
                        period = list(1995:2019),
                        clusters = list(~cluster),
                        strata = list(NULL),
                        id = list("unique_id"),
                        dob = list("wdob"),
                        intv = list("doi"),
                        weight = list("weight"),
                        varmethod = list("none"),
                        bhdata = mics_births_asfr,
                        bvars = list("cdob"),
                        counts = TRUE) %>%
    bind_rows %>%
    type.convert() %>%
    separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
    dplyr::filter(period <= survyear) %>%
    rename(value = asfr) %>%
    # rename(age_group = agegr) %>%
    mutate(survtype = "MICS",
           iso3 = iso3,
           variable = "asfr"
    ) %>%
    left_join(get_age_groups() %>% select(age_group, age_group_label), by=c("agegr" = "age_group_label")) %>%
    select(-agegr)

  mics_asfr_plot_nat <- Map(calc_asfr, mics_wm_asfr,
                            by = list(~survey_id),
                            tips = list(c(0,15)),
                            agegr= list(3:10*5),
                            period = list(1995:2019),
                            clusters = list(~cluster),
                            strata = list(NULL),
                            id = list("unique_id"),
                            dob = list("wdob"),
                            intv = list("doi"),
                            weight = list("weight"),
                            varmethod = list("none"),
                            bhdata = mics_births_asfr,
                            bvars = list("cdob"),
                            counts = TRUE) %>%
    dplyr::bind_rows() %>%
    utils::type.convert() %>%
    tidyr::separate(col = survey_id, into = c(NA, "survyear", NA), sep = c(3,7), remove = FALSE, convert = TRUE) %>%
    dplyr::filter(period <= survyear) %>%
    dplyr::rename(value = asfr) %>%
    # rename(age_group = agegr) %>%
    dplyr::mutate(survtype = "MICS",
           iso3 = iso3,
           area_id = iso3,
           variable = "asfr"
    ) %>%
    dplyr::left_join(get_age_groups() %>% select(age_group, age_group_label), by=c("agegr" = "age_group_label")) %>%
    dplyr::select(-agegr)


  mics_tfr_plot <- mics_asfr_plot %>%
    dplyr::group_by(survey_id, period, area_id) %>%
    dplyr::summarise(value = 5*sum(value)) %>%
    dplyr::mutate(variable = "tfr") %>%
    dplyr::filter(value > 0.5)

  mics_tfr_plot_nat <- mics_asfr_plot_nat %>%
    dplyr::group_by(survey_id, period, area_id) %>%
    dplyr::summarise(value = 5*sum(value)) %>%
    dplyr::mutate(variable = "tfr") %>%
    dplyr::filter(value > 0.5)

  missing_data <- mics_asfr %>%
    dplyr::group_by(survey_id, tips) %>%
    dplyr::summarise(births = sum(births)) %>%
    dplyr::filter(births < 5) %>%
    dplyr::group_by(survey_id) %>%
    dplyr::group_split()

  mics_plot <- dplyr::bind_rows(mics_tfr_plot, mics_tfr_plot_nat, mics_asfr_plot, mics_asfr_plot_nat)

  i <- 0

  while(i < length(missing_data)) {
    i <- i+1

    mics_asfr <- mics_asfr %>%
      dplyr::filter(!(survey_id == unique(missing_data[[i]]$survey_id) & tips %in% unique(missing_data[[i]]$tips)))

    years <- unique(dplyr::filter(mics_asfr, survey_id == unique(missing_data[[i]]$survey_id))$period)

    mics_plot <- mics_plot %>%
      dplyr::filter(!(survey_id == unique(missing_data[[i]]$survey_id) & !period %in% years))

  }

  out <- list()
  out$mics_asfr <- mics_asfr
  out$mics_plot <- mics_plot

  out
}

#' Calculate fertility rates from DHS, MIS, AIS data
#' @export
calculate_dhs_fertility <- function(iso3, surveys, clusters, areas_wide) {

  mc.cores <- if(.Platform$OS.type == "windows") 1 else parallel::detectCores()

  cluster_list <- clusters %>%
    dplyr::rename(area_id = geoloc_area_id) %>%
    dplyr::group_by(survey_id) %>%
    dplyr::group_split()

  names(cluster_list) <- surveys$survey_id

  ir <- get_fertility_surveys(surveys)
  names(ir$ir) <- names(cluster_list)
  names(ir$nrow_ir) <- names(cluster_list)

  if(iso3 == 'ETH') {

    adjust_eth_months <- function(ir) {

      col_positions <- grep("v011|v008|^b3\\_[0-9]*", colnames(ir))
      ir <- ir %>%
        dplyr::mutate(dplyr::across(col_positions, ~{.x+92}))
    }

    ir$ir <- lapply(ir$ir, adjust_eth_months)

  }

  dat <- map_ir_to_areas(ir$ir, cluster_list)

  cluster_level <- clusters %>%
    dplyr::filter(!is.na(geoloc_area_id)) %>%
    tidyr::separate(geoloc_area_id, into = c(NA, "area_level", NA), sep = c(4,5), remove = FALSE, convert = TRUE) %>%
    dplyr::mutate(area_level = ifelse(geoloc_area_id == iso3, 0, area_level)) %>%
    dplyr::group_by(area_level) %>%
    dplyr::group_split()

  get_admin1_clusters <- function(cluster_level) {
    lvl <- unique(cluster_level$area_level)

    if(lvl > 1) {
      cluster_level %>%
        dplyr::left_join(areas_wide %>%
                    sf::st_drop_geometry() %>%
                    dplyr::select(area_id1, paste0("area_id", lvl)) %>%
                    dplyr::distinct(),
                  by = c("geoloc_area_id" = paste0("area_id", lvl))) %>%
        dplyr::rename(area_id = area_id1) %>%
        dplyr::select(survey_id, cluster_id, area_id)
    } else {
      cluster_level %>%
        dplyr::rename(area_id = geoloc_area_id) %>%
        dplyr::select(survey_id, cluster_id, area_id)
    }
  }

  cluster_list_admin1 <- lapply(cluster_level, get_admin1_clusters) %>%
    dplyr::bind_rows() %>%
    dplyr::ungroup() %>%
    dplyr::group_by(survey_id) %>%
    dplyr::group_split()

  names(cluster_list_admin1) <- surveys$survey_id

  dat_admin1 <- map_ir_to_areas(ir$ir, cluster_list_admin1, single_tips = FALSE)
  dat_admin1$ir <- lapply(dat_admin1$ir, zap_labels)

  message("Calculating district-level ASFR")

  asfr <- Map(calc_asfr, dat$ir,
              by = list(~survey_id + survtype + survyear + area_id),
              tips = list(c(0:15)),
              agegr= list(3:10*5),
              period = list(1995:2020),
              strata = list(NULL),
              varmethod = list("none"),
              counts = TRUE) %>%
    dplyr::bind_rows() %>%
    utils::type.convert() %>%
    dplyr::filter(period <= survyear) %>%
    # rename(age_group = agegr) %>%
    dplyr::mutate(iso3 = iso3) %>%
    dplyr::left_join(get_age_groups() %>% dplyr::select(age_group, age_group_label), by = c("agegr" = "age_group_label")) %>%
    dplyr::select(-agegr)

  message("Calculating aggregate fertility rates")

  asfr_plot <- Map(calc_asfr, dat_admin1$ir,
                   by = list(~survey_id + survtype + survyear + area_id),
                   tips = dat_admin1$tips_surv,
                   agegr= list(3:10*5),
                   period = list(1995:2020),
                   strata = list(NULL),
                   varmethod = list("none"),
                   counts = TRUE) %>%
    dplyr::bind_rows() %>%
    utils::type.convert() %>%
    dplyr::filter(period <= survyear) %>%
    dplyr::mutate(iso3 = iso3,
           variable = "asfr") %>%
    dplyr::rename(value = asfr) %>%
    dplyr::left_join(naomi::get_age_groups() %>% dplyr::select(age_group, age_group_label), by = c("agegr" = "age_group_label")) %>%
    dplyr::select(-agegr)

  asfr_plot_nat <- Map(calc_asfr, ir$ir,
                       # by = list(~survey_id + survtype + survyear + area_id),
                       tips = list(c(0,15)),
                       agegr= list(3:10*5),
                       period = list(1995:2020),
                       strata = list(NULL),
                       varmethod = list("none"),
                       counts = TRUE) %>%
    dplyr::bind_rows(.id = "survey_id") %>%
    tidyr::separate(survey_id, into = c(NA, "survyear", "survtype"), sep = c(3, 7), remove = FALSE, convert = TRUE) %>%
    utils::type.convert() %>%
    dplyr::filter(period<=survyear) %>%
    dplyr::mutate(iso3 = iso3,
           area_id = iso3,
           variable = "asfr") %>%
    dplyr::rename(value = asfr) %>%
    dplyr::left_join(naomi::get_age_groups() %>% dplyr::select(age_group, age_group_label), by = c("agegr" = "age_group_label")) %>%
    dplyr::select(-agegr)

  tfr_plot <- asfr_plot %>%
    dplyr::group_by(survey_id, period, area_id) %>%
    dplyr::summarise(value = 5 * sum(value)) %>%
    dplyr::mutate(variable = "tfr") %>%
    dplyr::filter(value > 0.5)

  tfr_plot_nat <- asfr_plot_nat %>%
    dplyr::group_by(survey_id, period, area_id) %>%
    dplyr::summarise(value = 5 * sum(value)) %>%
    dplyr::mutate(variable = "tfr") %>%
    dplyr::filter(value > 0.5)

  missing_data <- asfr %>%
    dplyr::group_by(survey_id, tips) %>%
    dplyr::summarise(births = sum(births)) %>%
    dplyr::filter(births < 5) %>%
    dplyr::group_by(survey_id) %>%
    dplyr::group_split()

  plot <- dplyr::bind_rows(tfr_plot, tfr_plot_nat, asfr_plot, asfr_plot_nat)

  i <- 0

  while(i < length(missing_data)) {
    i <- i+1

    asfr <- asfr %>%
      dplyr::filter(!(survey_id == unique(missing_data[[i]]$survey_id) & tips %in% unique(missing_data[[i]]$tips)))

    years <- unique(filter(asfr, survey_id == unique(missing_data[[i]]$survey_id))$period)

    plot <- plot %>%
      dplyr::filter(!(survey_id == unique(missing_data[[i]]$survey_id) & !period %in% years))


  }

  out <- list()
  out$asfr <- asfr
  out$plot <- plot
  out$nrow_ir <- ir$nrow_ir

  out

}
