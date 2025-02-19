#### 0 LIBRARIES AND PATHS ####################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)
library(INLA)
library(survey)
library(purrr)
library(raster)
library(sampling)
library(svymap)
inla.setOption("num.threads" = 8)
#### 0.1 File management ####
home_dir <- '~/'
if (!("Cluster" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
} else {
  #inla.setOption(enable.inla.argument.weights=TRUE)
}
setwd(paste0(home_dir, "Cluster/svyulm/"))
source("analysis/models.R")
out_res_dir <- "results/binary-simulations/"
dir.create(file.path(out_res_dir), showWarnings = FALSE)
sim_res_dir <- "analysis/binary-simulations/"
dir.create(file.path(sim_res_dir), showWarnings = FALSE)
sim_res_dir <- paste0(sim_res_dir, "results/")
dir.create(file.path(sim_res_dir), showWarnings = FALSE)
t0 <- Sys.time()
source("analysis/data_generating_functions.R")
source("analysis/sampling_functions.R")
#### 3 SIMULATION FUNCTIONS ####################################################
#### 3.1 Run single simulation ####
run_simulation <- function(formula, 
                           cov_dat, 
                           pop_config, 
                           sample_config,
                           n_cl_per_a, n_per_cl,
                           seed) {
  print(paste0(n_cl_per_a, ":",n_per_cl))
  set.seed(seed)
  pop_dat <- generate_binary_responses(cov_dat, config = pop_config)
  pop_dat$id_stra <- pop_dat$id_area
  pop_dat$domain <- as.character(pop_dat$id_area)
  pop_means <- pop_dat %>%
    group_by(domain) %>%
    summarize(pop_mean = mean(y)) %>%
    as.data.frame() 
  
  sample_dat <- generate_sample(sample_config, pop_dat, n_cl_per_a, n_per_cl)
  sample_dat <- sample_dat %>%
    group_by(id_stra) %>%
    mutate(wt_scl = wt / sum(wt) * n()) %>%
    ungroup()
  sample_dat$wt_placeholder = 1
  if (grepl("CLUS", sample_config)) {
    ids <- ~id_clus + id_unit
  } else {
    ids <- ~id_unit
  }
  
  unwt_des <- svydesign(ids = ids, strata = ~id_stra,
                        weights = ~wt_placeholder, data = sample_dat)
  wt_des <- svydesign(ids = ids, strata = ~id_stra,
                      weights = ~wt_scl, data = sample_dat)
  
  hajek_est <- get_direct(~y, ~domain, wt_des)
  greg_est <-  get_GREG(formula, ~domain, 
                        wt_des, pop_dat,
                        family = "quasibinomial") 
  Bayes_fit <-
    svymap::fit_psbayes(formula, 
                        ~domain, ~wt_placeholder,
                        svydes = unwt_des, 
                        family = "binomial",
                        rescale = F)
  Bayes_est<-
    svymap::get_psbayes_area_mean_binomial(formula, 
                                           pop_dat, ~domain, Bayes_fit$domain,
                                           Bayes_fit$draws) %>%
    mutate(method = "Bayes")

  psBayes_fit <-
    svymap::fit_psbayes(formula, 
                        ~domain, ~wt_scl, svydes = wt_des, 
                        family = "binomial",
                        rescale = T)
  
  psBayes_est<-
    svymap::get_psbayes_area_mean_binomial(formula, 
                                           pop_dat, ~domain, psBayes_fit$domain,
                                           psBayes_fit$draws, psBayes_fit$draws_rescale) 
  res <- bind_rows(hajek_est,
                   greg_est,
                   Bayes_est,
                   psBayes_est) 
  if (grepl("CLUS", sample_config)) {
    clus_dat <- sample_dat %>%
      group_by(id_clus) %>%
      mutate(y = sum(y),
             ntrials = n()) |>
      as.data.frame()
    clus_dat <- clus_dat[match(unique(clus_dat$id_clus), clus_dat$id_clus),]
    clus_wt_des <- svydesign(id = ~id_clus,
                             strata = ~id_stra, nest = T, 
                             weights = ~wt_scl, data = clus_dat)

    psBayes_bb_fit <-
      svymap::fit_psbayes(formula, 
                          ~domain, ~wt_scl, svydes = clus_wt_des, 
                          family = "betabinomial",  Ntrials = clus_dat$ntrials,
                          rescale = T)
    psBayes_bb_est<-
      svymap::get_psbayes_area_mean(formula, 
                                    pop_dat, ~domain, psBayes_bb_fit$domain,
                                    family = "betabinomial",
                                    psBayes_bb_fit$draws, 
                                    psBayes_bb_fit$draws_rescale) %>%
      mutate(method = paste0(method, "_bb"))
    res <- bind_rows(res, psBayes_bb_est)
    
  }
  

  res$id_sim = seed
  res <- left_join(res, pop_means, by = "domain")
  return(res)
}
#### 3.2 Example usage ####
# cov_dat <- generate_population_covariates(20, 50, 30)
# formula <- y ~ xo_unit + xo_clus
# res <- run_simulation(formula, cov_dat, "FEO", "SRS", n_cl_per_a = 3, n_per_cl = 3, seed = 1)
# ggplot(res, aes(x = pop_mean, y = est, color = method)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower, ymax = upper)) + 
#   geom_abline(slope = 1) + facet_wrap(~method)
#### 4 RUN SIMULATIONS #########################################################
for (i in 1:5) {
  cov_dat_path <- paste0("analysis/binary-simulations/binary-sims-cov-dat-", i, ".rds")
  if (!file.exists(cov_dat_path)) {
    cov_dat <- generate_population_covariates(20, 150, 30)
    saveRDS(cov_dat, cov_dat_path)
  } 
}

args = commandArgs(TRUE)
# supplied at the command line
seed = as.numeric(args[3])
pop_config <- as.character(args[1])
sample_config <- as.character(args[2])
sim_res_dir <- paste0(sim_res_dir, pop_config, "/")
dir.create(file.path(sim_res_dir), showWarnings = FALSE)
sim_res_dir <- paste0(sim_res_dir, sample_config, "/")
dir.create(file.path(sim_res_dir), showWarnings = FALSE)
res_list <- list()
# n_per_cl <- 5
# ncs <- c(6, 20)
# for (i in 1:1) {
#   cov_dat_path <- paste0("analysis/binary-simulations/binary-sims-cov-dat-", i, ".rds")
#   cov_dat <- readRDS(cov_dat_path)
#   for (nc in ncs) {
#     res <-  run_simulation(y ~ xo_unit,
#                            cov_dat, pop_config, sample_config,
#                            n_cl_per_a = nc, n_per_cl = n_per_cl, 
#                            seed = seed) %>%
#       mutate(pop_config = pop_config, 
#              sample_config = sample_config,
#              n_cl = nc,
#              n_per_cl = n_per_cl,
#              i = i) 
#     res_list <- c(res_list, list(res))
#   }
#   res <- bind_rows(res_list)
# }
# 
# res <- bind_rows(res_list)
# saveRDS(res, paste0(sim_res_dir, "res_", seed, ".rds"))

n_per_cl <- 5
ncs <- c(3)
for (i in 1:1) {
  cov_dat_path <- paste0("analysis/binary-simulations/binary-sims-cov-dat-", i, ".rds")
  cov_dat <- readRDS(cov_dat_path)
  for (nc in ncs) {
    res <-  run_simulation(y ~ xo_unit,
                           cov_dat, pop_config, sample_config,
                           n_cl_per_a = nc, n_per_cl = n_per_cl, 
                           seed = seed) %>%
      mutate(pop_config = pop_config, 
             sample_config = sample_config,
             n_cl = nc,
             n_per_cl = n_per_cl,
             i = i) 
    res_list <- c(res_list, list(res))
  }
  res <- bind_rows(res_list)
}

res <- bind_rows(res_list)
saveRDS(res, paste0(sim_res_dir, "small_sample_res_", seed, ".rds"))

# 
# pop_config <- "UC2"
# sample_config <- "CLUSINF"
# nc <- 20
# n_per_cl <- 5
# seed <- 4
# 
# res <-  run_simulation(y ~ xo_unit,
#                        cov_dat, pop_config, sample_config,
#                        n_cl_per_a = nc, n_per_cl = n_per_cl,
#                        seed = seed) %>%
#   mutate(pop_config = pop_config,
#          sample_config = sample_config,
#          n_cl = nc,
#          n_per_cl = n_per_cl)
# res %>% mutate(n_per_a = n_cl*n_per_cl) |> group_by(method, n_per_a) %>%
#   summarize(
#     rmse = sqrt(mean((est - pop_mean)^2)),
#     mae = mean(abs(est - pop_mean)),
#     mil = mean(upper - lower),
#     cov90 = mean(lower < pop_mean& upper > pop_mean)
#   )
# 
# 
# 
