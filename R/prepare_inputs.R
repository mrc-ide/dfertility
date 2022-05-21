#' Extract area, population, boundary files from target directory
#'
#' @param iso3_current A string of one or more ISO-3 codes
#' @param naomi_data_path A path to directory containing input data files
#' @return A named list of file paths

get_input_files <- function(iso3_current, naomi_data_path) {

  paths <- file.path(naomi_data_path, iso3_current, "data")

  files <- lapply(paths, function(paths) {

    files <- list.files(paths, full.names = TRUE)
    area <- files %>%
      stringr::str_subset(pattern = "areas.geojson") %>%
      stringr::str_subset(pattern = ".zip", negate = TRUE)

    if(population) {
      pop <- files %>%
        stringr::str_subset(pattern = "population")
      files <- c(area, pop)
      names(files) <- c("areas", "population")
    } else {
      files <- c(area)
      names(files) <- c("areas")
    }

    return(files)

  })

  names(files) <- iso3_current

  return(files)

}

#' Make area dataframes
#'
#' @description Create area hierarchy dataframes in the Naomi package format. The function will always return the area hierarchy in long format, with arguments additionally to return full and wide hierarchies
#' @param iso3_current A string of one or more ISO-3 codes
#' @param naomi_data_path A path to directory containing input data files
#' @param full Boolean to return the full area heirarchy
#' @param wide Boolean to return the area hierarchy in wide format
#' @export

get_areas <- function(iso3_current, naomi_data_path, full=FALSE, wide=TRUE) {

  files <- get_input_files(iso3_current, naomi_data_path)

  areas <- list()

  areas$areas_long <- lapply(files, "[[", "areas") %>%
    lapply(sf::read_sf) %>%
    lapply(function(x) {

      iso3_code <- x %>%
        dlpyr::filter(area_level == 0) %>%
        dlpyr::select(area_id) %>%
        unique() %>%
        .$area_id

      x <- x %>%
        dplyr::mutate(iso3 = iso3_code) %>%
        sf::st_drop_geometry() %>%
        dplyr::select(c("iso3", "area_id", "area_name", "area_level", "parent_area_id", "naomi_level"))

      return(x)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(iso3)

  if(full)
    areas$areas_full <- lapply(files, "[[", "areas") %>%
      lapply(st_read) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(iso3)

  if(wide)
    areas$areas_wide <- lapply(files, "[[", "areas") %>%
      lapply(sf::read_sf) %>%
      lapply(function(x) {naomi::spread_areas(as.data.frame(x))}) %>%
      lapply(function(x) {x %>% dplyr::mutate(iso3 = area_id0)}) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(iso3)

  return(areas)

}

#' Get populations
#' @description Returns a dataframe of population by district, five year age group, and sex
#' @param iso3_current A string of one or more ISO-3 codes
#' @param naomi_data_path A path to directory containing input data files
#' @export

get_populations <- function(iso3_current, naomi_data_path) {

  files <- get_input_files(iso3_current, naomi_data_path)

  area_population <- lapply(files, "[[", "population") %>%
    lapply(readr::read_csv) %>%
    lapply(dplyr::left_join, areas_long) %>%
    dlpyr::bind_rows() %>%
    dlpyr::mutate(period = naomi::year_labels(naomi:::calendar_quarter_to_quarter_id(calendar_quarter))) %>%
    dlpyr::select(iso3, "area_id" , "area_name", "source", "sex", "age_group", "population", "period") %>%
    dlpyr::arrange(iso3)
}

#' Get populations
#' @description Returns a dataframe of area boundaries
#' @param iso3_current A string of one or more ISO-3 codes
#' @param naomi_data_path A path to directory containing input data files
#' @export

get_boundaries <- function(iso3_current, naomi_data_path) {

  files <- get_input_files(iso3_current, naomi_data_path)

  area_boundaries <- lapply(files, "[[", "areas") %>%
    lapply(sf::read_sf) %>%
    lapply(function(x) {

      iso3_code <- x %>%
        dplyr::filter(area_level == 0) %>%
        dplyr::select(area_id) %>%
        unique() %>%
        .$area_id

      x <- x %>%
        dplyr::mutate(iso3 = iso3_code) %>%
        dplyr::select(-epp_level)

      return(x)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(iso3)
}
