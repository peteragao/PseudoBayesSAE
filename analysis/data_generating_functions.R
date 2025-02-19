#### 1 GENERATE POPULATION #####################################################
#### 1.1 Generate covariates ####
generate_population_covariates <- function(m, N_cl_per_area, N_per_cl) {
  N_cl <- m * N_cl_per_area
  N <- N_cl * N_per_cl
  cov_dat <- data.frame(
    id_area = rep(1:m, each = N / m),
    id_clus = rep(1:N_cl, each = N / N_cl),
    id_unit = 1:N
  ) %>%
    mutate(xo_unit = rnorm(N, mean = id_area / m, sd = 1),
           xo_clus = rep(rnorm(N_cl, mean = id_area / m, sd = 1),  
                         each = N / N_cl),
           xo_area = rep(rnorm(m, mean = id_area / m, sd = 1),  
                         each = N_cl_per_area * N_per_cl),
           size_unit = (id_area / m) + rexp(N, rate = 1/2)) %>%
    group_by(id_clus) %>%
    mutate(size_clus = sum(size_unit)) %>%
    ungroup() %>%
    group_by(id_area) %>%
    mutate(size_area = sum(size_unit)) %>%
    ungroup() %>%
    mutate(size_unit_scl = scale(size_unit)[,1],
           size_clus_scl = scale(size_clus)[,1],
           size_area_scl = scale(size_area)[,1])
  cov_dat 
}

#### 1.2 Generate random effects and responses ####
generate_binary_responses <- function(cov_dat,
                                      beta = c(0, 1, 0, .5, 1, 0),
                                      sd_area = 1,
                                      sd_clus = 1,
                                      config = c(NULL, "AO", "UU", "UC", 
                                                 "UA", "RA", "RC")[1]) {
  if (config == "AO") {
    beta <- c(0, 1, 0, 0, 0, 0, 0)
    sd_area <- 0
    sd_clus <- 0
  } else if (config == "UU") {
    beta <- c(0, 1, 0, 0, 2, 0, 0)
    sd_area <- 0
    sd_clus <- 0
  } else if (config == "UC") {
    beta <- c(0, 1, 0, 0, 0, 2, 0)
    sd_area <- 0
    sd_clus <- 0
  } else if (config == "UC2") {
      beta <- c(0, 1, 0, 0, 0, 1, 0)
      sd_area <- 0
      sd_clus <- .5
  } else if (config == "UA") {
    beta <- c(0, 1, 0, 0, 0, 0, 1)
    sd_area <- 0
    sd_clus <- 0
  } else if (config == "RA") {
    beta <- c(0, 1, 0, 0, 0, 0, 0)
    sd_area <- 1
    sd_clus <- 0
  } else if (config == "RC") {
    beta <- c(0, 1, 0, 0, 0, 0, 0)
    sd_area <- sqrt(.5)
    sd_clus <- sqrt(.5)
  }
  cov_mat <- 
    model.matrix(~ xo_unit + xo_clus + xo_area +
                   size_unit_scl + size_clus_scl + size_area_scl, 
                 cov_dat)
  pop_dat <- cov_dat %>%
    mutate(fe = (cov_mat %*% beta)[,1]) %>%
    group_by(id_area) %>%
    mutate(u = rep(rnorm(1, sd = sd_area), each = n()),
           ybar = mean(fe) + u) %>%
    ungroup() %>%
    group_by(id_clus) %>%
    mutate(v = rep(rnorm(1, sd = sd_clus), each = n())) %>%
    ungroup() %>%
    mutate(y = rbinom(n(), size = 1, prob = expit(fe + u + v)))
}
generate_gaussian_responses <- function(cov_dat,
                                        beta = c(0, 1, 0, .5, 1, 0),
                                        sd_area = 1,
                                        sd_clus = 1,
                                        sd_resp = 1,
                                        config = c(NULL, "AO", "UU", "UC", 
                                                   "UA", "RA", "RC")[1]) {
  if (config == "AO") {
    beta <- c(0, 1, 0, 0, 0, 0, 0)
    sd_area <- 0
    sd_clus <- 0
    sd_resp <- 1 
  } else if (config == "UU") {
    beta <- c(0, 1, 0, 0, 2, 0, 0)
    sd_area <- 0
    sd_clus <- 0
    sd_resp <- 1
  } else if (config == "UC") {
    beta <- c(0, 1, 0, 0, 0, 2, 0)
    sd_area <- 0
    sd_clus <- 0
    sd_resp <- 1
  } else if (config == "UC2") {
      beta <- c(0, 1, 0, 0, 0, 1, 0)
      sd_area <- 0
      sd_clus <- .5
 } else if (config == "UA") {
    beta <- c(0, 1, 0, 0, 0, 0, 1)
    sd_area <- 0
    sd_clus <- 0
    sd_resp <- 1
  } else if (config == "RA") {
    beta <- c(0, 1, 0, 0, 0, 0, 0)
    sd_area <- 1
    sd_clus <- 0
    sd_resp <- 1
  } else if (config == "RC") {
    beta <- c(0, 1, 0, 0, 0, 0, 0)
    sd_area <- sqrt(.5)
    sd_clus <- sqrt(.5)
    sd_resp <- 1
  }
  cov_mat <- 
    model.matrix(~ xo_unit + xo_clus + xo_area +
                   size_unit_scl + size_clus_scl + size_area_scl, 
                 cov_dat)
  pop_dat <- cov_dat %>%
    mutate(fe = (cov_mat %*% beta)[,1]) %>%
    group_by(id_area) %>%
    mutate(u = rep(rnorm(1, sd = sd_area), each = n()),
           ybar = mean(fe) + u) %>%
    ungroup() %>%
    group_by(id_clus) %>%
    mutate(v = rep(rnorm(1, sd = sd_clus), each = n())) %>%
    ungroup() %>%
    mutate(e = rnorm(n(), sd = sd_resp), 
           y = fe + u + v + e)
  
}
