convert_age_groups <- function(df) {

  df %>%
    dplyr::left_join(naomi::get_age_groups() %>%
                dplyr::select(age_group, age_group_label) %>%
                dplyr::rename(age_group_hold = age_group),
              by = c("age_group" = "age_group_label")) %>%
    dplyr::select(-age_group) %>%
    dplyr::rename(age_group = age_group_hold)

}

