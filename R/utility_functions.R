convert_age_groups <- function(df) {
  
  df %>%
    left_join(get_age_groups() %>%
                select(age_group, age_group_label) %>%
                rename(age_group_hold = age_group), 
              by=c("age_group" = "age_group_label")) %>%
    select(-age_group) %>%
    rename(age_group = age_group_hold)
  
}

