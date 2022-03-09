get_fertility_surveys <- function(surveys) {

  ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)

  ird <- ird %>%
    mutate(path = unlist(get_datasets(.))) %>%
    bind_rows()

  ir <- lapply(ird$path, readRDS) %>%
    lapply(function(x) {class(x) <- "data.frame"
    return(x)}) %>%
    Map(function(ir, surveys) {
      mutate(ir,
             surveyid = surveys$survey_id,
             country = surveys$CountryName,
             survyear = surveys$SurveyYear,
             survtype = surveys$SurveyType)
    }, ., group_split(surveys, SurveyId))

  ir


}

assign_tips <- function(ir, single_tips) {

  survey_type <- ir %>%
    lapply("[", "survtype") %>%
    lapply(unique) %>%
    bind_rows

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

  ir <- Map(ir_by_area, ir, cluster_list[names(ir)], n=1:length(ir), total=length(ir)) %>%
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
#' @export

assign_cluster_area <- function(clusters, areas_wide, area_level) {

  areas <- clusters %>%
    rename(area_id = geoloc_area_id) %>%
    left_join(areas_wide %>% select(area_id, paste0("area_id", area_level))) %>%
    select(-area_id) %>%
    rename(area_id = paste0("area_id", area_level))

  area_list <- areas %>%
    group_by(survey_id) %>%
    group_split(keep=TRUE)

  names(area_list) <- area_list %>%
    lapply("[", "survey_id") %>%
    lapply(unique) %>%
    bind_rows %>%
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
    filter(!is.na(area_id)) %>%
    group_split(area_id)

  return(ir_int)

}

#' Join individual recode survey datasets by cluster to area IDs
#' @param surveys Survey dataframe as created by \code{rdhs::dhs_surveys()}
#' @param cluster_areas Dataframe of clusters assigned to areas. Output of \code{\link{cluster_areas}}.
#' @param single_tips desc this
#' @return Nested list: \code{df$ir} and \code{df$tips_surv}. \code{df$ir} has one list item per survey/area combination and \code{df$tips_surv} has the corroborating TIPS value.
#' @export
#' @seealso \code{\link{cluster_areas}}

clusters_to_surveys <- function(iso3, surveys, cluster_areas, level, single_tips = TRUE) {

  ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)

  ird <- ird %>%
    mutate(path = unlist(get_datasets(.))) %>%
    bind_rows()

  ir <- lapply(ird$path, readRDS) %>%
    lapply(function(x) {class(x) <- "data.frame"
    return(x)}) %>%
    Map(function(ir, surveys) {
      mutate(ir,
             surveyid = surveys$survey_id,
             country = surveys$CountryName,
             survyear = surveys$SurveyYear,
             survtype = surveys$SurveyType)
    }, ., group_split(surveys, SurveyId))


  ## I think this mess was necessary because the order of cluster_areas was not always the same as the order for ir.. Double check. Otherwise the simple names(ir) <- names(cluster_areas) is sufficient
  # names(ir) <- ir %>%
  #   lapply("[", "surveyid") %>%
  #   lapply(unique) %>%
  #   bind_rows %>%
  #   left_join(clusters %>%
  #               select(survey_id, DHS_survey_id) %>%
  #               unique,
  #             by=c("surveyid" = "DHS_survey_id")) %>%
  #   select(-surveyid) %>%
  #   .$survey_id

  names(ir) <- names(cluster_areas)

  if(iso3 == "ETH") {

    cols_edit <- c("v008", "v011", "b3_01", "b3_02", "b3_03", "b3_04", "b3_05", "b3_06", "b3_07", "b3_08", "b3_09", "b3_10", "b3_11", "b3_12", "b3_13", "b3_14", "b3_15", "b3_16", "b3_17", "b3_18", "b3_19", "b3_20")

    ir <- ir %>%
      lapply(function(x) x %>% mutate_at(.vars = cols_edit, .funs = ~(.+92)))

  }


  if(level > 0) {

    ir <- Map(ir_by_area, ir, cluster_areas[names(ir)], n=1:length(ir), total=length(ir)) %>%
      unlist(recursive = FALSE)

    survey_type <- ir %>%
      lapply("[", "survtype") %>%
      lapply(unique) %>%
      bind_rows

    if(!single_tips)
      tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[survey_type$survtype]
    else
      tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[survey_type$survtype]

  } else {

    ir <- lapply(ir, function(x) {
      mutate(x, area_id = unique(surveys$iso3))
    })

    if(!single_tips)
      tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[surveys$SurveyType]
    else
      tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[surveys$SurveyType]

  }

  dat <- list()
  dat$ir <- ir
  dat$tips_surv <- tips_surv

  return(dat)

}

#' Calculate fertility rates from MICS data
#' @export
#'
calculate_mics_fertility <- function(iso3, mics_wm, mics_births_to_women) {

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

  #' For model:
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
    filter(period <= survyear) %>%
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
    filter(period <= survyear) %>%
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
    bind_rows %>%
    type.convert() %>%
    separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
    filter(period <= survyear) %>%
    rename(value = asfr) %>%
    # rename(age_group = agegr) %>%
    mutate(survtype = "MICS",
           iso3 = iso3,
           area_id = iso3,
           variable = "asfr"
    ) %>%
    left_join(get_age_groups() %>% select(age_group, age_group_label), by=c("agegr" = "age_group_label")) %>%
    select(-agegr)

  mics_tfr_plot <- Map(calc_tfr, mics_wm_asfr,
                       by = list(~area_id + survey_id),
                       tips = list(c(0,15)),
                       agegr= list(3:10*5),
                       period = list(1995:2020),
                       clusters = list(~cluster),
                       strata = list(NULL),
                       id = list("unique_id"),
                       dob = list("wdob"),
                       intv = list("doi"),
                       weight = list("weight"),
                       bhdata = mics_births_asfr,
                       bvars = list("cdob")) %>%
    bind_rows %>%
    type.convert() %>%
    separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
    filter(period <= survyear,
           tfr > 0.5) %>%
    rename(value = tfr) %>%
    # rename(age_group = agegr) %>%
    mutate(survtype = "MICS",
           iso3 = iso3,
           variable = "tfr"
    )

  mics_tfr_plot_nat <- Map(calc_tfr, mics_wm_asfr,
                           by = list(~survey_id),
                           tips = list(c(0,15)),
                           agegr= list(3:10*5),
                           period = list(1995:2020),
                           clusters = list(~cluster),
                           strata = list(NULL),
                           id = list("unique_id"),
                           dob = list("wdob"),
                           intv = list("doi"),
                           weight = list("weight"),
                           bhdata = mics_births_asfr,
                           bvars = list("cdob")) %>%
    bind_rows %>%
    type.convert() %>%
    separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
    filter(period <= survyear,
           tfr > 0.5) %>%
    rename(value = tfr) %>%
    # rename(age_group = agegr) %>%
    mutate(survtype = "MICS",
           iso3 = iso3,
           area_id = iso3,
           variable = "tfr"
    )

  missing_data <- mics_asfr %>%
    group_by(survey_id, tips) %>%
    summarise(births = sum(births)) %>%
    filter(births < 5) %>%
    group_by(survey_id) %>%
    group_split()

  mics_plot <- bind_rows(mics_tfr_plot, mics_tfr_plot_nat, mics_asfr_plot, mics_asfr_plot_nat)

  i <- 0

  while(i < length(missing_data)) {
    i <- i+1

    mics_asfr <- mics_asfr %>%
      filter(!(survey_id == unique(missing_data[[i]]$survey_id) & tips %in% unique(missing_data[[i]]$tips)))

    years <- unique(filter(mics_asfr, survey_id == unique(missing_data[[i]]$survey_id))$period)

    mics_plot <- mics_plot %>%
      filter(!(survey_id == unique(missing_data[[i]]$survey_id) & !period %in% years))

  }

  out <- list()
  out$mics_asfr <- mics_asfr
  out$mics_plot <- mics_plot

  out

}

#' Calculate fertility rates from DHS, MIS, AIS data
#' @export
#'
calculate_dhs_fertility <- function(iso3, surveys, clusters, areas_wide) {

  cluster_list <- clusters %>%
    rename(area_id = geoloc_area_id) %>%
    group_by(survey_id) %>%
    group_split

  names(cluster_list) <- surveys$survey_id

  ir <- get_fertility_surveys(surveys)
  names(ir) <- names(cluster_list)

  dat <- map_ir_to_areas(ir, cluster_list)

  cluster_list_admin1 <- clusters %>%
    left_join(areas_wide %>% st_drop_geometry, by=c("geoloc_area_id" = "area_id")) %>%
    rename(area_id = area_id1) %>%
    select(survey_id, cluster_id, area_id) %>%
    group_by(survey_id) %>%
    group_split

  names(cluster_list_admin1) <- surveys$survey_id

  dat_admin1 <- map_ir_to_areas(ir, cluster_list_admin1, single_tips = FALSE)
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
    bind_rows %>%
    type.convert %>%
    filter(period<=survyear) %>%
    # rename(age_group = agegr) %>%
    mutate(iso3 = iso3) %>%
    left_join(get_age_groups() %>% select(age_group, age_group_label), by=c("agegr" = "age_group_label")) %>%
    select(-agegr)

  message("Calculating aggregate fertility rates")

  asfr_plot <- Map(calc_asfr, dat_admin1$ir,
                   by = list(~survey_id + survtype + survyear + area_id),
                   tips = dat_admin1$tips_surv,
                   agegr= list(3:10*5),
                   period = list(1995:2020),
                   strata = list(NULL),
                   varmethod = list("none"),
                   counts = TRUE) %>%
    bind_rows %>%
    type.convert %>%
    filter(period<=survyear) %>%
    # rename(age_group = agegr) %>%
    mutate(iso3 = iso3,
           variable = "asfr") %>%
    rename(value = asfr) %>%
    left_join(get_age_groups() %>% select(age_group, age_group_label), by=c("agegr" = "age_group_label")) %>%
    select(-agegr)

  asfr_plot_nat <- Map(calc_asfr, ir,
                       # by = list(~survey_id + survtype + survyear + area_id),
                       tips = list(c(0,15)),
                       agegr= list(3:10*5),
                       period = list(1995:2020),
                       strata = list(NULL),
                       varmethod = list("none"),
                       counts = TRUE) %>%
    bind_rows(.id = "survey_id") %>%
    separate(survey_id, into=c(NA, "survyear", "survtype"), sep=c(3,7), remove=FALSE, convert=TRUE) %>%
    type.convert %>%
    filter(period<=survyear) %>%
    # rename(age_group = agegr) %>%
    mutate(iso3 = iso3,
           area_id = iso3,
           variable = "asfr") %>%
    rename(value = asfr) %>%
    left_join(get_age_groups() %>% select(age_group, age_group_label), by=c("agegr" = "age_group_label")) %>%
    select(-agegr)

  calc_tfr_wrapper <- function() {
    Map(calc_tfr, dat_admin1$ir,
        by = list(~survey_id + survtype + survyear + area_id),
        tips = dat_admin1$tips_surv,
        agegr= list(3:10*5),
        period = list(1995:2020),
        strata = list(NULL)) %>%
      bind_rows() %>%
      type.convert %>%
      filter(period<=survyear) %>%
      mutate(iso3 = iso3,
             variable = "tfr") %>%
      rename(value = tfr)
  }

  possibly_tfr <- possibly(calc_tfr_wrapper, otherwise = NULL)

  tfr_plot <- possibly_tfr()

  if(is.null(tfr_plot))
    tfr_plot <- asfr_plot %>%
      group_by(iso3, area_id, survey_id, survyear, survtype, period, age_group) %>%
      summarise(asfr = sum(births)/sum(pys)) %>%
      group_by(iso3, area_id, survey_id, survyear, survtype, period) %>%
      summarise(value = 5*sum(asfr)) %>%
      mutate(iso3 = iso3,
           variable = "tfr")

  # tfr_plot <- Map(calc_tfr, dat_admin1$ir,
  #                 by = list(~survey_id + survtype + survyear + area_id),
  #                 tips = dat_admin1$tips_surv,
  #                 agegr= list(3:10*5),
  #                 period = list(1995:2020),
  #                 strata = list(NULL)) %>%
  #   bind_rows() %>%
  #   type.convert %>%
  #   filter(period<=survyear) %>%
  #   mutate(iso3 = iso3,
  #          variable = "tfr") %>%
  #   rename(value = tfr)

  tfr_plot_nat <- Map(calc_tfr, ir,
                      # by = list(~survey_id + survtype + survyear + area_id),
                      tips = list(c(0,15)),
                      agegr= list(3:10*5),
                      period = list(1995:2020),
                      strata = list(NULL)) %>%
    bind_rows(.id = "survey_id") %>%
    separate(survey_id, into=c(NA, "survyear", "survtype"), sep=c(3,7), remove=FALSE, convert=TRUE) %>%
    type.convert %>%
    mutate(iso3 = iso3,
           variable = "tfr",
           area_id = iso3) %>%
    rename(value = tfr)

  missing_data <- asfr %>%
    group_by(survey_id, tips) %>%
    summarise(births = sum(births)) %>%
    filter(births < 5) %>%
    group_by(survey_id) %>%
    group_split()

  plot <- bind_rows(tfr_plot, tfr_plot_nat, asfr_plot, asfr_plot_nat)

  i <- 0

  while(i < length(missing_data)) {
    i <- i+1

    asfr <- asfr %>%
      filter(!(survey_id == unique(missing_data[[i]]$survey_id) & tips %in% unique(missing_data[[i]]$tips)))

    years <- unique(filter(asfr, survey_id == unique(missing_data[[i]]$survey_id))$period)

    plot <- plot %>%
      filter(!(survey_id == unique(missing_data[[i]]$survey_id) & !period %in% years))


  }

  out <- list()
  out$asfr <- asfr
  out$plot <- plot

  out

}
