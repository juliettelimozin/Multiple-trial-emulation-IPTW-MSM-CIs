CI_summary <- as.data.frame(cbind(0:4, t(pivot_coverage_ind[,,1,1] ))) %>% 
  mutate(Scenario = 1, Event_rate = 1) 
for (j in 1:3){
  for (i in 1:27){
    if(i == 1 & j == 1){ next }
    CI_summary_bis <- as.data.frame(cbind(0:4, t(pivot_coverage_ind[,,i,j] ))) %>% 
      mutate(Scenario = i, Event_rate = j) 
    CI_summary <- rbind(CI_summary, CI_summary_bis)
  }
}
colnames(CI_summary) <- c('Visit', 'Bootstrap', 'LEF outcome',
                          'LEF both','Sandwich','Jackknife MVN', 'Jackknife Wald', 'Scenario' , 'Event_rate')
CI_summary_long <- pivot_longer(CI_summary, cols = 2:7, names_to = 'CI_type', values_to = 'Coverage') %>% 
  dplyr::group_by(Visit, Scenario, Event_rate) %>% 
  dplyr::mutate(rank = dense_rank(desc(Coverage))) %>% 
  merge(check, by.x =  c('Event_rate', 'Scenario', 'Visit', 'CI_type'), by.y = c('j', 'i', 'Visit', 'CI_type')) %>% 
  dplyr::arrange(Event_rate, Scenario, Visit, rank)


CI_summary_old <- as.data.frame(cbind(0:4, t(pivot_coverage_ind_old[,,1,1] ))) %>% 
  mutate(Scenario = 1, Event_rate = 1) 
for (j in 1:3){
  for (i in 1:27){
    if(i == 1 & j == 1){ next }
    CI_summary_bis <- as.data.frame(cbind(0:4, t(pivot_coverage_ind_old[,,i,j] ))) %>% 
      mutate(Scenario = i, Event_rate = j) 
    CI_summary_old <- rbind(CI_summary_old, CI_summary_bis)
  }
}
colnames(CI_summary_old) <- c('Visit', 'Bootstrap', 'LEF outcome',
                          'LEF both','Sandwich','Jackknife MVN', 'Jackknife Wald', 'Scenario' , 'Event_rate')
CI_summary_long_old <- pivot_longer(CI_summary_old, cols = 2:7, names_to = 'CI_type', values_to = 'Coverage') %>% 
  dplyr::group_by(Visit, Scenario,Event_rate) %>% 
  dplyr::mutate(rank = dense_rank(desc(Coverage))) %>% 
  merge(check_old, by.x =  c('Event_rate', 'Scenario', 'Visit', 'CI_type'), by.y = c('j', 'i', 'Visit', 'CI_type')) %>% 
  dplyr::arrange(Event_rate, Scenario, Visit, rank)


CI_summary_long$old_rank <- CI_summary_long_old$CI_type
CI_summary_long$coverage_old <- CI_summary_long_old$Coverage
CI_summary_long$mean_old <- CI_summary_long_old$mean
CI_summary_long$q1_old <- CI_summary_long_old$q1
CI_summary_long$q3_old <- CI_summary_long_old$q3


save(CI_summary_long, file = 'CI_summary_long.rda')
write.csv(CI_summary_long,"CI_summary_long_same_seed.csv", row.names = FALSE)

CI_summary_by_type <- CI_summary_long %>% 
  dplyr::select(Event_rate, Scenario, Visit, CI_type, Coverage, rank) %>% 
  merge(CI_summary_long_old %>% select(-rank), by = c('Event_rate', 'Scenario', 'Visit', 'CI_type')) %>% 
  merge(check, by.x =  c('Event_rate', 'Scenario', 'Visit', 'CI_type'), by.y = c('j', 'i', 'Visit', 'CI_type')) %>% 
  merge(check_old,by.x =  c('Event_rate', 'Scenario', 'Visit', 'CI_type'), by.y = c('j', 'i', 'Visit', 'CI_type'), all.x = TRUE) %>% 
  dplyr::arrange(Event_rate, Scenario, Visit, rank) %>% 
  dplyr::select( -rank)
  

names(CI_summary_by_type) <- c('Event_rate', 'Scenario', 'Visit', 'CI_type', 'Coverage_fixed', 'Coverage_error', 'mean_se_fixed',
                               'q1_se_fixed', 'q3_se_fixed', 'mean_se_old', 'q1_se_old', 'q3_se_old')


write.csv(CI_summary_by_type,"CI_summary_by_type_same_seed.csv", row.names = FALSE)




