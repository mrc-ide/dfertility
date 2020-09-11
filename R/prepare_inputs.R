#' Extract area, population, boundary files from target directory
#'
#' @param iso3_current A string of one or more ISO-3 codes
#' @param naomi_data_path A path to directory containing input data files
#' @return A named list of file paths

get_input_files <- function(iso3_current, naomi_data_path) {

  paths <- file.path(naomi_data_path, iso3_current, "data")

  files <- lapply(paths, function(paths) {

    files <- list.files(paths, full.names = TRUE)
    area <- files %>% str_subset(pattern = "areas.geojson") %>% str_subset(pattern = ".zip", negate=TRUE)

    if(population) {
      pop <- files %>% str_subset(pattern = "population")
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
#'
#' @param iso3_current A string of one or more ISO-3 codes
#' @param naomi_data_path A path to directory containing input data files
#' @param full Boolean to return the full area heirarchy
#' @param wide Boolean to return the area hierarchy in wide format
#' @export

get_areas <- function(iso3_current, naomi_data_path, full=FALSE, wide=TRUE) {

  files <- get_input_files(iso3_current, naomi_data_path)

  areas <- list()

  areas$areas_long <- lapply(files, "[[", "areas") %>%
    lapply(read_sf) %>%
    lapply(function(x) {

      iso3_code <- x %>%
        filter(area_level == 0) %>%
        select(area_id) %>%
        unique %>%
        .$area_id

      x <- x %>%
        mutate(iso3 = iso3_code) %>%
        st_drop_geometry() %>%
        select(c("iso3", "area_id", "area_name", "area_level", "parent_area_id", "naomi_level"))

      return(x)
    }) %>%
    bind_rows %>%
    arrange(iso3)

  if(full)
    areas$areas_full <- lapply(files, "[[", "areas") %>%
      lapply(st_read) %>%
      bind_rows %>%
      arrange(iso3)

  if(wide)
    areas$areas_wide <- lapply(files, "[[", "areas") %>%
      lapply(read_sf) %>%
      lapply(function(x) {spread_areas(as.data.frame(x))}) %>%
      lapply(function(x) {x %>% mutate(iso3 = area_id0)}) %>%
      bind_rows %>%
      arrange(iso3)

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
    lapply(read_csv) %>%
    lapply(left_join, areas_long) %>%
    bind_rows %>%
    mutate(period = year_labels(naomi:::calendar_quarter_to_quarter_id(calendar_quarter))) %>%
    select(iso3, "area_id" , "area_name", "source", "sex", "age_group", "population", "period") %>%
    arrange(iso3)
}

#' Get populations
#' @description Returns a dataframe of area boundaries
#' @param iso3_current A string of one or more ISO-3 codes
#' @param naomi_data_path A path to directory containing input data files
#' @export

get_boundaries <- function(iso3_current, naomi_data_path) {

  files <- get_input_files(iso3_current, naomi_data_path)

  area_boundaries <- lapply(files, "[[", "areas") %>%
    lapply(read_sf) %>%
    lapply(function(x) {

      iso3_code <- x %>%
        filter(area_level == 0) %>%
        select(area_id) %>%
        unique %>%
        .$area_id

      x <- x %>%
        mutate(iso3 = iso3_code) %>%
        select(-epp_level)

      return(x)
    }) %>%
    bind_rows %>%
    arrange(iso3)
}
