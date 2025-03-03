CI_summary <- as.data.frame(cbind(0:4, t(pivot_coverage_ind[,,1,1] ))) %>% 
  mutate(Scenario = 1) 

for (i in 2:9){
  CI_summary_bis <- as.data.frame(cbind(0:4, t(pivot_coverage_ind[,,i,1] ))) %>% 
    mutate(Scenario = i) 
  CI_summary <- rbind(CI_summary, CI_summary_bis)
}
colnames(CI_summary) <- c('Visit', 'Bootstrap', 'LEF_outcome',
                          'LEF_both','Sandwich','Jackknife_MVN', 'Jackknife_Wald', 'Scenario' )
CI_summary_long <- pivot_longer(CI_summary, cols = Bootstrap:Jackknife_Wald, names_to = 'CI_type', values_to = 'Coverage') %>% 
  dplyr::group_by(Visit, Scenario) %>% 
  dplyr::mutate(rank = dense_rank(desc(Coverage))) %>% 
  dplyr::arrange(Scenario, Visit, rank)


CI_summary_old <- as.data.frame(cbind(0:4, t(pivot_coverage_ind_old[,,1,1] ))) %>% 
  mutate(Scenario = 1) 

for (i in 2:9){
  CI_summary_bis <- as.data.frame(cbind(0:4, t(pivot_coverage_ind_old[,,i,1] ))) %>% 
    mutate(Scenario = i) 
  CI_summary_old <- rbind(CI_summary_old, CI_summary_bis)
}
colnames(CI_summary_old) <- c('Visit', 'Bootstrap', 'LEF_outcome',
                          'LEF_both','Sandwich','Jackknife_MVN', 'Jackknife_Wald', 'Scenario' )
CI_summary_long_old <- pivot_longer(CI_summary_old, cols = Bootstrap:Jackknife_Wald, names_to = 'CI_type', values_to = 'Coverage') %>% 
  dplyr::group_by(Visit, Scenario) %>% 
  dplyr::mutate(rank = dense_rank(desc(Coverage))) %>% 
  dplyr::arrange(Scenario, Visit, rank)

CI_summary_long$old_rank <- CI_summary_long_old$CI_type
CI_summary_long$coverage_old <- CI_summary_long_old$Coverage

save(CI_summary_long, file = 'CI_summary_long.rda')
write.csv(CI_summary_long,"CI_summary_long.csv", row.names = FALSE)


