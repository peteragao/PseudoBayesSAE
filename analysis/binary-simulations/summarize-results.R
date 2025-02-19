library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)
library(INLA)
library(survey)
library(purrr)
library(rgdal)
library(raster)
library(rstan)
library(sampling)
library(ggrepel)
library(kableExtra)
#### FILE MANAGEMENT ####
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}
setwd(paste0(home_dir, "Cluster/svyulm/"))
sample_configs <- c(
  "SRS",
  "SSINF",
  "CLUSINF"
)
pop_configs <- c(
  "UC2"
)
smry_list <- list()
for (pop_config in pop_configs) {
  for (sample_config in sample_configs) {
    current_sim_res_dir <- paste0("analysis/binary-simulations/results/",
                                  pop_config, "/", sample_config, "/")
    sim_res_list <- list()
    ctr = 0
    m <- 0
    for (k in 1:1000) {
      if (file.exists( 
        paste0(current_sim_res_dir,
               "res_", k, ".rds")
      )) {
        if (read <- T) {
          res <- readRDS(paste0(current_sim_res_dir,
                                "res_", k, ".rds")) %>%
            mutate(n_per_a = n_cl*n_per_cl)
          sim_res_list <- c(sim_res_list, list(res))
        }
        ctr <- ctr + 1
      } else {
        m <- k
        
      }
    }
    print(paste0("check ", m))
    print(paste0(pop_config, "-", sample_config, ": missing ", 1000-ctr))
    if (ctr > 0) {
      res <- bind_rows(sim_res_list)
      
      smry <- res %>% group_by(method, n_per_a, i) %>%
        summarize(
          rmse = sqrt(mean((est - pop_mean)^2)),
          mae = mean(abs(est - pop_mean)),
          mil = mean(upper - lower),
          cov90 = mean(lower < pop_mean& upper > pop_mean)
        ) %>%
        ungroup() %>%
        mutate(pop_config = pop_config,
               sample_config = sample_config)
      smry_list <- c(smry_list, list(smry))
      ggplot(smry, aes(x = n_per_a, y = rmse, color = method)) + geom_line()
      ggsave(paste0("analysis/binary-simulations/results/figures/",
                    pop_config, "-", sample_config, "_rmse.pdf"), width = 8, height = 7)
      ggplot(smry, aes(x = n_per_a, y = cov90, color = method)) + geom_line() +
      geom_hline(yintercept = .9) 
      ggsave(paste0("analysis/binary-simulations/results/figures/",
                    pop_config, "-", sample_config, "_cov90.pdf"), width = 8, height = 7)
    }
  }
}
all_smry <- bind_rows(smry_list)
ggplot(all_smry, aes(x = n_per_a, y = rmse, color = method)) +
  geom_line(size = 1) + 
  facet_wrap(pop_config~sample_config, scales = "free") + 
  theme_bw() 

ggplot(all_smry , aes(x = n_per_a, y = cov90, color = method)) +
  geom_line(size = 1) +
  geom_hline(yintercept = .9) + 
  facet_wrap(pop_config~sample_config, scales = "free") + 
  theme_bw() 


all_smry <- all_smry %>%
  filter(method %in% c("HAJEK", "GREG",
                       "Bayes", "psBayes", "psBayes_rscl",
                       "psBayes_bb", "psBayes_rscl_bb")[1:5]) %>%
  mutate(method = factor(method, 
                         levels = c("HAJEK", "GREG",
                                    "Bayes", "psBayes", "psBayes_rscl",
                                    "psBayes_bb", "psBayes_rscl_bb"))) %>%
  mutate(sample_config = factor(sample_config,
                                levels = c("SRS", "SSINF", "CLUSINF")))
levels(all_smry$sample_config) <- list("SRS"  = "SRS", 
                                       "PPS1" = "SSINF",
                                       "PPS2" = "CLUSINF")
levels(all_smry$method) <- list("Hájek"  = "HAJEK", 
                                "GREG" = "GREG",
                                "Unwt" = "Bayes",
                                "Wt" = "psBayes",
                                "WtRscl" = "psBayes_rscl",
                                "WtBB" = "psBayes_bb",
                                "WtRsclBB" = "psBayes_rscl_bb")

fmt_tbl <- function(res, label, caption, file = NULL, methods = unique(res$method), 
                    rows = NULL, bold_min = F) {
  res <- res %>%
    #   filter(method %in% methods) %>%
    #   mutate(method = 
    #            pub_names$publication[match(method, pub_names$internal)]) %>%
    #   arrange(match(method, pub_order)) %>%
    dplyr::select(sample_config, method, rmse,mae, mil, cov90) %>%
    mutate(rmse = rmse * 100, mae = mae * 100, mil = mil * 100, cov90 = cov90 * 100) %>%
    setNames(c("Design",
               "Method", 
               paste0("RMSE (x 100)"),
               "MAE (x 100)",
               "MIL (x 100)",
               "90% Int. Cov.")) %>%
    knitr::kable(digits = c(0, 0, 1, 1, 1, 0), format = "latex", booktabs = T,
                 linesep = "", label = label, caption = caption) %>%
    kableExtra::collapse_rows(columns = 1) %>%
    kableExtra::kable_styling(font_size = 7)
  if (!is.null(rows)) {
    #res <- res %>%  kableExtra::row_spec(rows, hline_after = T) 
  }
  if (!is.null(file)) {
    res %>% writeLines(file)
  }
  res
}
res <- all_smry %>% filter(n_per_a == 30, i == 1) %>% arrange(pop_config, sample_config, method)
capt <- paste0("Averaged evaluation metrics of estimators of area level means",
               " across 1,000 binary response simulations for SRS, PPS1, and", 
               " PPS2 designs for a sample size of thirty units per area.")
lbl <- "bin-sim-smry"
fmt_tbl(res, label = lbl, caption = capt, 
        file =  'paper/figures/binary_simulation_summary.tex')
lbl <- "bin-sim-smry-lg"
capt <- paste0("Averaged evaluation metrics of estimators of area level means",
               " across 1,000 binary response simulations for SRS, PPS1, and", 
               " PPS2 designs for a sample size of one hundred units per area.")
res <- all_smry %>% filter(n_per_a == 100, i == 1) %>% arrange(pop_config, sample_config, method)
fmt_tbl(res, label = lbl, caption = capt, 
        'paper/figures/binary_simulation_summary_large_n.tex')

#### Summaries of model fit ####
r2 <- numeric()
for (k in 1:30) {
  pop_config <- "UC2"
  i <- 1
  cov_dat_path <- paste0("analysis/binary-simulations/binary-sims-cov-dat-", i, ".rds")
  cov_dat <- readRDS(cov_dat_path)
  
  source("analysis/data_generating_functions.R")
  source("analysis/sampling_functions.R")
  pop_dat <- generate_binary_responses(cov_dat, config = pop_config)
  
  
  #### Model fit - Cohen's R2 ####
  pop_fit <- glm(y ~ xo_unit, data = pop_dat, family = "quasibinomial")
  r2[k] <- with(summary(pop_fit), 1 - deviance/null.deviance)
}



pop_fit <- glm(y ~ xo_unit + size_clus_scl, data = pop_dat, family = "quasibinomial")
with(summary(pop_fit), 1 - deviance/null.deviance)

smry_list <- list()
for (pop_config in pop_configs) {
  for (sample_config in sample_configs) {
    current_sim_res_dir <- paste0("analysis/binary-simulations/results/",
                                  pop_config, "/", sample_config, "/")
    sim_res_list <- list()
    ctr = 0
    m <- 0
    for (k in 1:1000) {
      if (file.exists( 
        paste0(current_sim_res_dir,
               "small_sample_res_", k, ".rds")
      )) {
        if (read <- T) {
          res <- readRDS(paste0(current_sim_res_dir,
                                "small_sample_res_", k, ".rds")) %>%
            mutate(n_per_a = n_cl*n_per_cl)
          sim_res_list <- c(sim_res_list, list(res))
        }
        ctr <- ctr + 1
      } else {
        m <- k
        
      }
    }
    print(paste0("check ", m))
    print(paste0(pop_config, "-", sample_config, ": missing ", 1000-ctr))
    if (ctr > 0) {
      res <- bind_rows(sim_res_list)
      
      smry <- res %>% group_by(method, n_per_a, i) %>%
        summarize(
          rmse = sqrt(mean((est - pop_mean)^2)),
          mae = mean(abs(est - pop_mean)),
          mil = mean(upper - lower),
          cov90 = mean(lower < pop_mean& upper > pop_mean)
        ) %>%
        ungroup() %>%
        mutate(pop_config = pop_config,
               sample_config = sample_config)
      smry_list <- c(smry_list, list(smry))
     
    }
  }
}
all_smry <- bind_rows(smry_list)
ggplot(all_smry, aes(x = n_per_a, y = rmse, color = method)) +
  geom_line(size = 1) + 
  facet_wrap(pop_config~sample_config, scales = "free") + 
  theme_bw() 

ggplot(all_smry , aes(x = n_per_a, y = cov90, color = method)) +
  geom_line(size = 1) +
  geom_hline(yintercept = .9) + 
  facet_wrap(pop_config~sample_config, scales = "free") + 
  theme_bw() 


all_smry <- all_smry %>%
  filter(method %in% c("HAJEK", "GREG",
                       "Bayes", "psBayes", "psBayes_rscl",
                       "psBayes_bb", "psBayes_rscl_bb")[1:5]) %>%
  mutate(method = factor(method, 
                         levels = c("HAJEK", "GREG",
                                    "Bayes", "psBayes", "psBayes_rscl",
                                    "psBayes_bb", "psBayes_rscl_bb"))) %>%
  mutate(sample_config = factor(sample_config,
                                levels = c("SRS", "SSINF", "CLUSINF")))
levels(all_smry$sample_config) <- list("SRS"  = "SRS", 
                                       "PPS1" = "SSINF",
                                       "PPS2" = "CLUSINF")
levels(all_smry$method) <- list("Hájek"  = "HAJEK", 
                                "GREG" = "GREG",
                                "Unwt" = "Bayes",
                                "Wt" = "psBayes",
                                "WtRscl" = "psBayes_rscl",
                                "WtBB" = "psBayes_bb",
                                "WtRsclBB" = "psBayes_rscl_bb")

fmt_tbl <- function(res, label, caption, file = NULL, methods = unique(res$method), 
                    rows = NULL, bold_min = F) {
  res <- res %>%
    #   filter(method %in% methods) %>%
    #   mutate(method = 
    #            pub_names$publication[match(method, pub_names$internal)]) %>%
    #   arrange(match(method, pub_order)) %>%
    dplyr::select(sample_config, method, rmse,mae, mil, cov90) %>%
    mutate(rmse = rmse * 100, mae = mae * 100, mil = mil * 100, cov90 = cov90 * 100) %>%
    setNames(c("Design",
               "Method", 
               paste0("RMSE (x 100)"),
               "MAE (x 100)",
               "MIL (x 100)",
               "90% Int. Cov.")) %>%
    knitr::kable(digits = c(0, 0, 1, 1, 1, 0), format = "latex", booktabs = T,
                 linesep = "", label = label, caption = caption) %>%
    kableExtra::collapse_rows(columns = 1) %>%
    kableExtra::kable_styling(font_size = 7)
  if (!is.null(rows)) {
    #res <- res %>%  kableExtra::row_spec(rows, hline_after = T) 
  }
  if (!is.null(file)) {
    res %>% writeLines(file)
  }
  res
}
lbl <- "bin-sim-smry-sml"
capt <- paste0("Averaged evaluation metrics of estimators of area level means",
               " across 1,000 binary response simulations for SRS, PPS1, and", 
               " PPS2 designs for a sample size of fifteen units per area.")
res <- all_smry %>% filter(n_per_a == 15, i == 1) %>% arrange(pop_config, sample_config, method)
fmt_tbl(res, label = lbl, caption = capt, 
        'paper/figures/binary_simulation_summary_small_n.tex')

















