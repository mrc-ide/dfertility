#' Download MICS surveys
#' @param iso3 ISO3 code
#' @export

create_surveys_mics <- function(iso3, mics_indicators) {

  sharepoint <- spud::sharepoint$new(Sys.getenv("SHAREPOINT_URL"))
  folder <- sharepoint$folder(site = Sys.getenv("SHAREPOINT_SITE"), path = Sys.getenv("MICS_ORDERLY_PATH"))

  mics_file_names <- folder$list() %>%
    dplyr::filter(stringr::str_detect(name, tolower(iso3))) %>%
    .$name %>%
    sort()

  mics_survey_names <- toupper(stringr::str_replace(mics_file_names, ".rds", ""))

  mics_file_names <- lapply(tolower(unique(dplyr::filter(mics_indicators, survey_id %in% mics_survey_names)$survey_id)),
         grep, x = mics_file_names, value = TRUE
  ) %>%
    unlist() %>%
    sort()


  rename_datasets_key <- mics_indicators %>%
    dplyr::filter(survey_id %in% mics_survey_names) %>%
    dplyr::filter(label == "dataset name") %>%
    dplyr::arrange(survey_id, filetype) %>%
    dplyr::group_by(survey_id) %>%
    dplyr::group_split()

  paths <- file.path("sites", Sys.getenv("SHAREPOINT_SITE"), Sys.getenv("MICS_ORDERLY_PATH"), mics_file_names)

  files <- lapply(paths, spud::sharepoint_download, sharepoint_url = Sys.getenv("SHAREPOINT_URL"))

  mics_dat <- lapply(files, readRDS)

  extract_rename_datasets <- function(mics_dat, rename_datasets_key) {

    mics_dat <- mics_dat[rename_datasets_key$value]

    names(mics_dat) <- rename_datasets_key$id

    return(mics_dat)

  }

  mics_dat <- Map(extract_rename_datasets, mics_dat, rename_datasets_key)

  names(mics_dat) <- mics_indicators %>%
    dplyr::filter(survey_id %in% mics_survey_names) %>%
    dplyr::distinct(survey_id) %>%
    .$survey_id %>%
    sort()

  return(mics_dat)

}

#' Filter MICS datasets
#' @description Filter MICS household, women, and birth history datasets to key variables, and rename to ensure consistent column names between surveys
#' @export

filter_mics <- function(dat, mics_indicators, survey_id_i) {

  indicators <- dplyr::filter(mics_indicators,
                       survey_id == survey_id_i,
                       label != "dataset name"
  )

  wm <- dat$wm
  colnames(wm) <- tolower(colnames(wm))
  wm <- wm %>%
    dplyr::select(dplyr::filter(indicators, filetype == "wm")$value)
  colnames(wm) <- dplyr::filter(indicators, filetype == "wm")$id

  bh <- dat$bh
  colnames(bh) <- tolower(colnames(bh))
  bh <- bh %>%
    dplyr::select(dplyr::filter(indicators, filetype == "bh")$value)
  colnames(bh) <- dplyr::filter(indicators, filetype == "bh")$id


  hh <- dat$hh
  colnames(hh) <- tolower(colnames(hh))
  hh <- hh %>%
    dplyr::select(dplyr::filter(indicators, filetype == "hh")$value)
  colnames(hh) <- dplyr::filter(indicators, filetype == "hh")$id

  df <- list()
  df$wm <- wm %>%
    dplyr::mutate(survey_id = survey_id_i) %>%
    dplyr::filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
    dplyr::arrange(cluster, hh_number, line_number) %>%
    dplyr::group_by(cluster, hh_number, line_number) %>%
    dplyr::mutate(unique_id = dplyr::cur_group_id())

  df$bh <- bh %>%
    dplyr::mutate(survey_id = survey_id_i)

  df$hh <- hh %>%
    dplyr::mutate(survey_id = survey_id_i,
           area_level = as.numeric(dplyr::filter(indicators, id == "mics_area_level")$value)
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
    dplyr::bind_rows(.id = "survey_id")

  hh <- mics_dat %>%
    lapply("[[", "hh") %>%
    lapply(function(x) {
      x %>%
        dplyr::left_join(data.frame(mics_area_name = attr(x$mics_area_name, "labels"),
                             mics_area_name_label = stringr::str_to_title(
                               names(attr(x$mics_area_name, "labels")))
                             )
                  ) %>%
        dplyr::select(-mics_area_name) %>%
        dplyr::mutate(mics_area_name_label = stringr::str_remove_all(mics_area_name_label, "\\u0093|\\u0094"),
               mics_area_name_label = stringr::str_trim(mics_area_name_label))

    }) %>%
    dplyr::bind_rows(.id = "survey_id")


  bh <- mics_dat %>%
    lapply("[[", "bh") %>%
    dplyr::bind_rows(.id = "survey_id")

  df <- list()
  df$wm <- wm
  df$hh <- hh
  df$bh <- bh

  return(df)
}

#' Join household datasets with area hierarchy
#' @export

join_survey_areas <- function(fertility_mics_data, areas, warn=FALSE) {

  errfun <- if (warn)
    warning
  else stop

  areas <- sf::st_drop_geometry(areas)

  dat <- fertility_mics_data$hh

  lvl <- dat %>%
    dplyr::distinct(survey_id, area_level)

  area_survey_level <- areas %>%
    dplyr::left_join(lvl)

  dat_merge <- dat %>%
    dplyr::full_join(area_survey_level %>%
                dplyr::select(survey_id, area_id, area_name, area_level) %>%
                dplyr::filter(!is.na(survey_id)),
              by = c("mics_area_name_label" = "area_name", "area_level", "survey_id")
    )


  if (any(is.na(dat_merge$area_id))) {

    missing_areas <- dat_merge %>%
      dplyr::filter(is.na(area_id)) %>%
      dplyr::select(survey_id, mics_area_name_label, area_id) %>%
      dplyr::distinct() %>%
      dplyr::rename(mics_area_name = mics_area_name_label)

    errfun("\n\nSurvey regions were not matched to areas:\n",
         paste0(utils::capture.output(missing_areas), collapse = "\n"),
         "\n\nThis must be corrected \n \n"
    )
  }

  if(any(is.na(dat_merge$cluster))) {

    missing_areas <- dat_merge %>%
      dplyr::filter(is.na(cluster)) %>%
      dplyr::select(survey_id, mics_area_name_label, area_id) %>%
      dplyr::distinct() %>%
      dplyr::rename(area_name = mics_area_name_label)

    warning("\n\nAreas were not found in MICS survey:\n",
            paste0(utils::capture.output(missing_areas), collapse = "\n"),
            "\n\nThis may be because the survey did not sample these regions"
          )

  }

  fertility_mics_data$hh <- dat %>%
    dplyr::left_join(area_survey_level %>%
                dplyr::select(survey_id, area_id, area_name, area_level),
              by = c("mics_area_name_label" = "area_name", "area_level", "survey_id")
    )

  fertility_mics_data

}

#' Transform survey datasets into inputs for calc_asfr
#' @export
make_asfr_inputs <- function(mics_survey_areas, mics_survey_data) {

  dat <- mics_survey_areas
  df <- list()

  df$wm <- dat$wm %>%
    dplyr::left_join(dat$hh %>% dplyr::select(survey_id, cluster, hh_number, area_id))

  df$births_to_women <- dat$wm %>%
    dplyr::select(survey_id, cluster, hh_number, line_number, unique_id) %>%
    dplyr::left_join(dat$bh %>% dplyr::select(survey_id, cluster, hh_number, line_number, cdob)) %>%
    dplyr::select(survey_id, unique_id, cdob) %>%
    dplyr::filter(!is.na(cdob))

  df
}
