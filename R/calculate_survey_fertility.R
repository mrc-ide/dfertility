#' Aggregate district level cluster dataframe
#' @description Aggregate a dataframe of survey clusters geolocated to an administrative to any higher administrative level
#' @param clusters Dataframe of geolocated survey clusters
#' @param areas_wide Area hierarchy in wide format
#' @param area_level Integer referring to administrative level, where 0 is national, 1 is provincial etc.
#' @return Dataframe list of clusters by administrative area, one list item per survey
#' @seealso \code{\link{get_areas}}
#' @export

assign_cluster_area <- function(clusters, areas_wide, area_level) {

  iso3_current <- unique(clusters$iso3)

  if (area_level == "naomi") {

    int <- areas_long %>%
      filter(iso3 == iso3_current)

    area_level <- unique(int$area_level[int$naomi_level == TRUE])

  }

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

#' Join individual recode survey datasets by cluster to area IDs
#' @param surveys
#' @param cluster_areas Dataframe of clusters assigned to areas. Output of \code{\link{cluster_areas}}.
#' @param single_tips <<describe>>
#' @return Nested list: \code{df$ir} and \code{df$tips_surv}. \code{df$ir} has one list item per survey/area combination and \code{df$tips_surv} has the corroborating TIPS value.
#' @export
#' @seealso \code{\link{cluster_areas}}

clusters_to_surveys <- function(surveys, cluster_areas, single_tips = TRUE) {

  level <- areas_long$area_level[areas_long$area_id == cluster_areas[[1]]$area_id[[1]]]

  ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)

  ird <- ird %>%
    mutate(path = unlist(get_datasets(.))) %>%
    bind_rows()

  ir <- lapply(ird$path, readRDS) %>%
    lapply(function(x) {class(x) <- "data.frame"
    return(x)}) %>%
    Map(function(ir, surveys) {
      mutate(ir,
             surveyid = surveys$SurveyId,
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

  if(unique(surveys$iso3) == "ETH") {

    cols_edit <- c("v008", "v011", "b3_01", "b3_02", "b3_03", "b3_04", "b3_05", "b3_06", "b3_07", "b3_08", "b3_09", "b3_10", "b3_11", "b3_12", "b3_13", "b3_14", "b3_15", "b3_16", "b3_17", "b3_18", "b3_19", "b3_20")

    ir <- ir %>%
      lapply(function(x) x %>% mutate_at(.vars = cols_edit, .funs = ~(.+92)))

  }


  if(level > 0) {

    ir <- Map(ir_by_area2, ir, cluster_areas[names(ir)], n=1:length(ir), total=length(ir)) %>%
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