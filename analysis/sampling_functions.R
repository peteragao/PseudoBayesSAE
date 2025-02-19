#### GENERATE SAMPLE #####################################################
#### 1.1. Stratified simple random sampling ####
generate_sample_SRS <- function(pop_dat, n_cl_per_a, n_per_cl) {
  n_per_a <- n_cl_per_a * n_per_cl
  pop_dat <- pop_dat %>% 
    group_by(id_stra) %>%
    mutate(pi_2 = n_per_a / n(),
           wt_2 = 1 / pi_2,
           s_2 = srswor(n_per_a, n())) %>%
    ungroup() %>%
    mutate(wt = wt_2, wt_1 = 1)
  sample_dat <- pop_dat %>%
    filter(s_2 == 1)
  return(sample_dat)
}
#### 1.2. Stratified single stage uninformative ####
generate_sample_SS <- function(pop_dat, n_cl_per_a, n_per_cl) {
  n_per_a <- n_cl_per_a * n_per_cl
  pop_dat <- pop_dat %>% 
    group_by(id_stra) %>%
    mutate(su = rexp(n(), rate = 1 / 2),
           pi_2 = inclusionprobabilities(su - min(su) + 1, n_per_a),
           wt_2 = 1 / pi_2,
           s_2 = UPmaxentropy(pi_2)) %>%
    ungroup() %>%
    mutate(wt = wt_2, wt_1 = 1)
  sample_dat <- pop_dat %>%
    filter(s_2 == 1)
  return(sample_dat)
}

#### 1.3. Stratified single stage informative ####
generate_sample_SSINF <- function(pop_dat, n_cl_per_a, n_per_cl) {
  n_per_a <- n_cl_per_a * n_per_cl
  pop_dat <- pop_dat %>% 
    group_by(id_stra) %>%
    mutate(pi_2 = inclusionprobabilities(size_unit_scl - min(size_unit_scl) + 1, n_per_a),
           wt_2 = 1 / pi_2,
           s_2 = UPmaxentropy(pi_2)) %>%
    ungroup() %>%
    mutate(wt = wt_2, wt_1 = 1)
  sample_dat <- pop_dat %>%
    filter(s_2 == 1)
  return(sample_dat)
}
#### 1.4. Stratified two stage cluster ####
generate_sample_CLUS <- function(pop_dat, n_cl_per_a, n_per_cl) {
  pop_dat <- pop_dat %>%
    group_by(id_clus) %>%  
    mutate(pi_2 = n_per_cl / n(),
           wt_2 = 1 / pi_2,
           s_2 = srswor(n_per_cl, n())) %>%
    ungroup() %>%
    group_by(id_stra, id_clus) %>%
    nest() %>%
    ungroup() %>% 
    group_by(id_stra) %>%
    mutate(pi_1 = n_cl_per_a / n(),
           wt_1 = 1 / pi_1,
           s_1 = srswor(n_cl_per_a, n())) %>%
    ungroup() %>%
    unnest(data) %>%
    ungroup() %>%
    mutate(wt = wt_1, wt_2 = 1)
  sample_dat <- pop_dat %>%
    filter(s_1 == 1 & s_2 == 1)
  return(sample_dat)
}

#### 1.5. Stratified two stage cluster informative ####
generate_sample_CLUSINF <- function(pop_dat, n_cl_per_a, n_per_cl) {
  pop_dat <- pop_dat %>%
    mutate(su = size_unit_scl - min(size_unit_scl) + 1,
           sc = size_clus_scl - min(size_clus_scl) + 1) %>%
    group_by(id_clus) %>%  
    mutate(pi_2 = inclusionprobabilities(su, n_per_cl),
           wt_2 = 1 / pi_2,
           s_2 = UPmaxentropy(pi_2)) %>%
    ungroup() %>%
    group_by(id_stra, id_clus, sc) %>%
    nest() %>%
    ungroup() %>% 
    group_by(id_stra) %>%
    mutate(pi_1 = inclusionprobabilities(sc, n_cl_per_a),
           wt_1 = 1 / pi_1,
           s_1 = UPmaxentropy(pi_1)) %>%
    ungroup() %>%
    unnest(data) %>%
    mutate(wt = wt_1 * wt_2) %>%
    dplyr::select(-c(su, sc))
  sample_dat <- pop_dat %>%
    filter(s_1 == 1 & s_2 == 1)
  return(sample_dat)
}
generate_sample <- function(sample_config, pop_dat, n_cl_per_a, n_per_cl) {
  if (sample_config == "SS") {
    return(generate_sample_SS(pop_dat, n_cl_per_a, n_per_cl))
  } else if (sample_config == "SSINF") {
    return(generate_sample_SSINF(pop_dat, n_cl_per_a, n_per_cl))
  } else if (sample_config == "CLUS") {
    return(generate_sample_CLUS(pop_dat, n_cl_per_a, n_per_cl))
  } else if (sample_config == "CLUSINF") {
    return(generate_sample_CLUSINF(pop_dat, n_cl_per_a, n_per_cl))
  }  else if (sample_config == "SRS") {
    return(generate_sample_SRS(pop_dat, n_cl_per_a, n_per_cl))
  }
}
#### 2 Summary visualizations of samples ####

#### 3 Example usage ####
# pop_dat$id_stra <- pop_dat$id_area
# sample_dat_SRS <- generate_sample_SRS(pop_dat, 2, 5)
# sample_dat_SS <- generate_sample_SS(pop_dat, 2, 5)
# sample_dat_SSINF <- generate_sample_SSINF(pop_dat, 2, 5)
# sample_dat_CLUS <- generate_sample_CLUS(pop_dat, 2, 5)
# sample_dat_CLUSINF <- generate_sample_CLUSINF(pop_dat, 2, 5)
# 
# # informative sampling
# mean(sample_dat_SRS$y)
# mean(sample_dat_SS$y)
# mean(sample_dat_SSINF$y)
# mean(sample_dat_CLUS$y)
# mean(sample_dat_CLUSINF$y)
# 
# # design effects
# sd(sample_dat_SRS$y)
# sd(sample_dat_SS$y)
# sd(sample_dat_SSINF$y)
# sd(sample_dat_CLUS$y)
# sd(sample_dat_CLUSINF$y)
