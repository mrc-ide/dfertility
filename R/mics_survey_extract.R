#' Download MICS surveys
#' @param iso3 ISO3 code
#' @export

create_surveys_mics <- function(iso3) {

  sharepoint <- spud::sharepoint$new(Sys.getenv("SHAREPOINT_URL"))
  folder <- sharepoint$folder(site = Sys.getenv("SHAREPOINT_SITE"), path = Sys.getenv("MICS_ORDERLY_PATH"))

  iso3_mics <- folder$list() %>%
    filter(str_detect(name, tolower(iso3))) %>%
    .$name

  mics_survey_names <- str_replace(iso3_mics, ".rds", "")

  paths <- file.path("sites", Sys.getenv("SHAREPOINT_SITE"), Sys.getenv("MICS_ORDERLY_PATH"), iso3_mics)

  files <- lapply(paths, sharepoint_download, sharepoint_url = Sys.getenv("SHAREPOINT_URL"))

  mics_dat <- lapply(files, readRDS) %>%
    lapply("[", c("wm", "bh", "hh"))

  names(mics_dat) <- toupper(mics_survey_names)

  return(mics_dat)

}

#' Filter MICS datasets
#' @description Filter MICS household, women, and birth history datasets to key variables, and rename to ensure consistent column names between surveys
#' @export

filter_mics <- function(dat, mics_indicators, survey_id_i) {

  indicators <- filter(mics_indicators, survey_id == survey_id_i)

  wm <- dat$wm
  colnames(wm) <- tolower(colnames(wm))
  wm <- wm %>%
    select(filter(indicators, filetype == "wm")$value)
  colnames(wm) <- filter(indicators, filetype == "wm")$id

  bh <- dat$bh
  colnames(bh) <- tolower(colnames(bh))
  bh <- bh %>%
    select(filter(indicators, filetype == "bh")$value)
  colnames(bh) <- filter(indicators, filetype == "bh")$id


  hh <- dat$hh
  colnames(hh) <- tolower(colnames(hh))
  hh <- hh %>%
    select(filter(indicators, filetype == "hh")$value)
  colnames(hh) <- filter(indicators, filetype == "hh")$id

  df <- list()
  df$wm <- wm %>%
    mutate(survey_id = survey_id_i) %>%
    filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
    arrange(cluster, hh_number, line_number) %>%
    mutate(unique_id = group_indices(., cluster, hh_number, line_number))

  df$bh <- bh %>%
    mutate(survey_id = survey_id_i)

  df$hh <- hh %>%
    mutate(survey_id = survey_id_i,
           area_level = as.numeric(filter(indicators, id == "mics_area_level")$value)
    )

  return(df)


}

#' Transform MICS dataframes
#' @description Convert lists by survey to lists by dataset type
#' @export
transform_mics <- function(mics_survey_data, mics_indicators) {

  mics_dat <- Map(filter_mics, mics_survey_data, list(mics_indicators), names(mics_survey_data))

  wm <- mics_dat %>%
    lapply("[[", "wm") %>%
    bind_rows(.id = "survey_id")

  hh <- mics_dat %>%
    lapply("[[", "hh") %>%
    lapply(function(x) {
      x %>%
        left_join(data.frame(mics_area_name = attr(x$mics_area_name, "labels"), mics_area_name_label = str_to_title(names(attr(x$mics_area_name, "labels")))))

    }) %>%
    bind_rows(.id = "survey_id")

  bh <- mics_dat %>%
    lapply("[[", "bh") %>%
    bind_rows(.id = "survey_id")

  df <- list()
  df$wm <- wm
  df$hh <- hh
  df$bh <- bh

  return(df)

}

#' Join household datasets with area hierarchy
#' @export

join_survey_areas <- function(fertility_mics_data, areas) {

  areas <- areas %>% st_drop_geometry()

  dat <- fertility_mics_data$hh

  lvl <- unique(dat$area_level)

  dat_merge <- dat %>%
    full_join(areas %>%
                filter(area_level ==lvl) %>%
                select(area_id, area_name, area_level),
              by=c("mics_area_name_label" = "area_name", "area_level")
    )


  if (any(is.na(dat_merge$area_id))) {

    missing_areas <- dat_merge %>%
      filter(is.na(area_id)) %>%
      select(survey_id, mics_area_name_label, area_id) %>%
      distinct %>%
      rename(mics_area_name = mics_area_name_label)

    stop("Survey regions were not matched to areas:\n",
         paste0(capture.output(missing_areas), collapse = "\n"))
  }

  if(any(is.na(dat_merge$mics_area_name_label))) {

    missing_areas <- dat_merge %>%
      filter(is.na(area_id)) %>%
      select(survey_id, mics_area_name_label, area_id) %>%
      distinct %>%
      rename(mics_area_name = mics_area_name_label)

    warning("Areas did not have matching survey regions:\n",
            paste0(capture.output(missing_areas), collapse = "\n"))

  }

  fertility_mics_data$hh <- dat %>%
    left_join(areas %>%
                filter(area_level ==lvl) %>%
                select(area_id, area_name, area_level),
              by=c("mics_area_name_label" = "area_name", "area_level")
    )

  fertility_mics_data

}

#' Transform survey datasets into inputs for calc_asfr
#' @export
make_asfr_inputs <- function(mics_survey_areas, mics_survey_data) {

  dat <- mics_survey_areas
  df <- list()

  df$wm <- dat$wm %>%
    left_join(dat$hh %>% select(survey_id, cluster, hh_number, area_id))

  df$births_to_women <- dat$wm %>%
    select(survey_id, cluster, hh_number, line_number, unique_id) %>%
    left_join(dat$bh %>% select(survey_id, cluster, hh_number, line_number, cdob)) %>%
    select(survey_id, unique_id, cdob) %>%
    filter(!is.na(cdob))

  df

}
