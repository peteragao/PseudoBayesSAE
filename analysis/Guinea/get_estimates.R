library(survey)
library(svymap)
library(INLA)
library(ggplot2)
library(sf)
library(terra)
library(spdep)
library(dplyr)
library(tidyr)
library(haven)
gadm_abbrev <- "GIN"
pop_abbrev <- "gin"
country <- "Guinea"
#### 1 GADM DATA #################################################################
poly_layer_adm0 <- paste('gadm41', gadm_abbrev,
                         '0', sep = "_")
poly_layer_adm1 <- paste('gadm41', gadm_abbrev,
                         '1', sep = "_")
poly_layer_adm2 <- paste('gadm41', gadm_abbrev,
                         '2', sep = "_")
poly_path <- "data/Guinea/gadm41_GIN_shp/"

poly_adm0 <- st_read(dsn = poly_path,
                     layer = as.character(poly_layer_adm0)) 
# use encoding to read special characters
poly_adm1 <-  st_read(dsn = poly_path,
                      layer = as.character(poly_layer_adm1)) 
if(sum(grepl(paste('gadm41', gadm_abbrev,
                   '2', sep = "_"), list.files(poly_path))) != 0){
  poly_adm2 <-  st_read(dsn = poly_path,
                        layer = as.character(poly_layer_adm2)) 
}
if (exists("poly_adm2")) {
  st_crs(poly_adm0) <- st_crs(poly_adm1)  <- st_crs(poly_adm2)
}else {
  st_crs(poly_adm0) <- st_crs(poly_adm1)
}
poly_adm2$NAME_1 <- as.character(poly_adm2$NAME_1)
poly_adm1$NAME_1 <- as.character(poly_adm1$NAME_1)
#### 1.1 Create adjacency matrices and tables ####
if(exists("poly_adm1")){
  admin1_mat <- poly2nb(poly_adm1)
  admin1_mat <- nb2mat(admin1_mat, zero.policy = TRUE)
  colnames(admin1_mat) <- 
    rownames(admin1_mat) <-
    paste0("admin1_", 1:dim(admin1_mat)[1])
  admin1_names <- data.frame(GADM = poly_adm1$NAME_1,
                             Internal = rownames(admin1_mat))
}else{
  message("There is no Admin1 polygon file.")
}
if(exists("poly_adm2")){
  admin2_mat <- poly2nb(poly_adm2)
  admin2_mat <- nb2mat(admin2_mat, zero.policy = TRUE)
  colnames(admin2_mat) <- 
    rownames(admin2_mat) <- 
    paste0("admin2_", 1:dim(admin2_mat)[1])
  admin2_names <- data.frame(GADM = poly_adm2$NAME_2,
                             Internal = rownames(admin2_mat))
}else{
  message("There is no Admin2 polygon file.")
}

#### 2 DHS DATA ################################################################
svy_dat <- readRDS("data/Guinea/Guinea_DHS_2018_MCV.rds")
svy_dat <- svy_dat %>%
  group_by(admin1_char) %>%
  mutate(wt_scl = wt / sum(wt) * n()) %>%
  ungroup()
# access raster
access <- rast("data/Nigeria/2015_accessibility_to_cities_v1.0.tif")

# nighttime lights raster
ntl <- rast(paste0("data/", country, "/Covariates/",
                   pop_abbrev,"_viirs_100m_2016.tif"))

svy_dat_coords <- st_coordinates(svy_dat)
svy_dat$access <- terra::extract(access, vect(svy_dat_coords))[, 2]

svy_dat$access <- ifelse(is.na(svy_dat$access), mean(svy_dat$access, na.rm = T), svy_dat$access)
svy_dat$l1a <- log(1 + svy_dat$access)
svy_dat$ntl <- terra::extract(ntl, vect(svy_dat_coords))[, 2]
svy_dat$ntl <- ifelse(is.na(svy_dat$ntl), mean(svy_dat$ntl, na.rm = T), svy_dat$ntl)
svy_dat$urban <- as.factor(svy_dat$urban)
svy_dat$wt_placeholder = 1

#### 3 POPULATION DATA #########################################################
pop_dat <- readRDS("data/Guinea/pix_tbl_1km/gin_pix_tbl.rds")
pop_dat$admin2_char <-
  admin2_names$Internal[match(pop_dat$admin2_name, admin2_names$GADM)]
pop_dat$urban <- factor(ifelse(pop_dat$urban, "U", "R"))




#### 4 GET ESTIMATES ###########################################################
#### 4.1 Estimation functions ####
get_direct<- function(response, domain, sample_des, CI = .90) {
  ests <- svyby(response, domain, sample_des, svymean)
  out <- data.frame(
    domain = ests[, 1],
    est = ests[, 2],
    median = ests[, 2],
    var = ests[, 3] ^ 2
  ) 
  out$lower <- qnorm((1-CI)/2) * sqrt(out$var) + out$est
  out$upper <- qnorm(1 - (1-CI)/2) * sqrt(out$var) + out$est
  out$method <- "HAJEK"
  return(out)
}
# temp <- get_direct(~y, ~id_area, sample_des)

get_GREG <- function(formula, domain, sample_des, newdata, 
                     family = "gaussian",
                     CI = .90, newweights = NULL) {
  
  domain_name <- as.character(domain)[2]
  sample_des$variables$domain_internal <- sample_des$variables[[domain_name]]
  newdata$domain_internal <- newdata[[domain_name]]
  
  working_fit <- svyglm(update(formula, .~. + domain_internal),
                        family = family, sample_des)
  pop_unit_ests <- as.vector(predict(working_fit, newdata,
                                     type = "response")) 
  if (is.null(newweights)) {
    newweights = rep(1, nrow(newdata))
  }
  pred_nums <-
    aggregate(pop_unit_ests * newweights, 
              list(domain_internal = as.vector(newdata$domain_internal)),
              sum)
  
  pred_denoms <-
    aggregate(newweights, 
              list(domain_internal = as.vector(newdata$domain_internal)),
              sum)
  area_ests <- data.frame(domain_internal = pred_nums[, 1], 
                          est = pred_nums[, 2] / pred_denoms[, 2])
  y <- model.response(model.frame(formula, sample_des$variables))
  
  sample_des$variables$res <- 
    y -
    as.vector(predict(working_fit, sample_des$variables, type = "response")) 
  res_ht <- svyby(~res, ~domain_internal, sample_des, svymean)
  out_dat <- left_join(area_ests, res_ht, by = "domain_internal")
  
  
  out_dat$domain <- out_dat$domain_internal
  
  out_dat$est = out_dat[, 2] + out_dat$res
  out_dat$median <- out_dat$est
  out_dat$var = out_dat$se ^ 2
  out_dat$method = "GREG"
  out_dat$lower = out_dat$est + qnorm((1-CI)/2) * out_dat$se
  out_dat$upper = out_dat$est + qnorm(1 - (1-CI)/2) * out_dat$se
  out_dat <- dplyr::select(out_dat, domain, est, median, var, lower, upper, method)
  return(arrange(out_dat, domain))
}

#### 4.2 Fit models ####
clus_dat <- svy_dat %>%
  group_by(cluster) %>%
  mutate(mcv_yes = sum(mcv_yes),
         ntrials = n()) 
clus_dat <- clus_dat[match(unique(clus_dat$cluster), clus_dat$cluster),]
clus_wt_des <- svydesign(id = ~cluster,
                      strata = ~stratum, nest=T, 
                      weights = ~wt_scl, data=clus_dat)
unwt_des <- svydesign(id = ~cluster + hshold,
                      strata = ~stratum, nest=T, 
                      weights = ~wt_placeholder, data=svy_dat)
wt_des <- svydesign(id = ~cluster + hshold,
                    strata = ~stratum, nest=T, 
                    weights = ~wt_scl, data=svy_dat)

hajek_est <- get_direct(~mcv_yes, ~admin2_char, wt_des)
greg_est <-  get_GREG(mcv_yes ~  urban + ntl + l1a,
                      ~admin2_char, wt_des, pop_dat,
                      family = "quasibinomial", 
                      newweights = pop_dat$svy_pop) 

system.time(
  psBayes_bb_fit <-
    fit_psbayes(mcv_yes ~ urban + ntl + l1a, 
                ~admin2_char, ~wt_placeholder, svydes = clus_wt_des, 
                family = "betabinomial", Ntrials = clus_dat$ntrials,
                rescale = T)  
)
system.time(
  Bayes_fit <-
    fit_psbayes(mcv_yes ~ urban + ntl + l1a, 
                ~admin2_char, ~wt_placeholder, svydes = unwt_des, 
                family = "binomial",
                rescale = F)  
)
system.time(
  psBayes_fit <-
    fit_psbayes(mcv_yes ~ urban + ntl + l1a, 
                ~admin2_char, ~wt_scl, svydes = wt_des, 
                family = "binomial",
                rescale = T)
  
)
system.time(
  Bayes_sp_fit <- 
    fit_psbayes(mcv_yes ~ urban + ntl + l1a, 
                ~admin2_char, ~wt_placeholder, 
                adj_mat = admin2_mat,
                svydes = unwt_des, 
                family = "binomial",
                rescale = F)
)
system.time(
  psBayes_sp_fit <- 
    fit_psbayes(mcv_yes ~ urban + ntl + l1a, 
                ~admin2_char, ~wt_scl,
                adj_mat = admin2_mat,
                svydes = wt_des, 
                family = "binomial",
                rescale = T)
)

Bayes_est <- get_psbayes_area_mean_binomial(
  mcv_yes ~ urban + ntl + l1a, 
  pop_dat, ~admin2_char, Bayes_fit$domain,
  Bayes_fit$draws,
  newweights = pop_dat$svy_pop
) %>%
  mutate(method = paste0("Bayes"))
Bayes_sp_est <- get_psbayes_area_mean_binomial(
  mcv_yes ~ urban + ntl + l1a, 
  pop_dat, ~admin2_char, Bayes_sp_fit$domain,
  Bayes_sp_fit$draws,
  newweights = pop_dat$svy_pop
)  %>%
  mutate(method = stringr::str_replace(method, "psBayes", "Bayes")) %>% 
  mutate(method = paste0(method, "_sp"))
psBayes_est <- get_psbayes_area_mean_binomial(
  mcv_yes ~ urban + ntl + l1a, 
  pop_dat, ~admin2_char, psBayes_fit$domain,
  psBayes_fit$draws, psBayes_fit$draws_rescale,
  newweights = pop_dat$svy_pop
) 

psBayes_bb_est <- get_psbayes_area_mean(
  mcv_yes ~ urban + ntl + l1a, 
  pop_dat, ~admin2_char, psBayes_bb_fit$domain,
  family = "betabinomial",
  psBayes_bb_fit$draws, psBayes_bb_fit$draws_rescale,
  newweights = pop_dat$svy_pop
) %>%
  mutate(method = paste0(method, "_bb"))
psBayes_sp_est <- get_psbayes_area_mean_binomial(
  mcv_yes ~ urban + ntl + l1a, 
  pop_dat, ~admin2_char, psBayes_sp_fit$domain,
  psBayes_sp_fit$draws, psBayes_sp_fit$draws_rescale,
  newweights = pop_dat$svy_pop
)  %>%
  mutate(method = paste0(method, "_sp"))

all_est <- bind_rows(hajek_est,
                     greg_est, 
                     Bayes_est,
                     Bayes_sp_est,
                     psBayes_est,
                     psBayes_sp_est) %>%
  left_join(admin2_names, by = c("domain" = "Internal"))

saveRDS(all_est, "results/Guinea/adm2_mcv_ests.rds")
#### 5 FIGURES #################################################################
#### 5.1 MAPS OF EA LOCATIONS ##################################################
all_est <- readRDS("results/Guinea/adm2_mcv_ests.rds")
all_est[all_est$GADM == "Kouroussa", 2:6] <- NA

all_est <- all_est %>%
  filter(method %in% c("HAJEK", "GREG",
                       "Bayes", "psBayes", "psBayes_rscl",
                       "psBayes_bb", "psBayes_rscl_bb")[1:5]) %>%
  mutate(method = factor(method, 
                         levels = c("HAJEK", "GREG",
                                    "Bayes", "psBayes", "psBayes_rscl",
                                    "Bayes_sp", "psBayes_sp", "psBayes_rscl_sp",
                                    "psBayes_bb", "psBayes_rscl_bb")))
levels(all_est$method) <- list("Hájek"  = "HAJEK", 
                               "GREG" = "GREG",
                               "Unwt" = "Bayes",
                               "Wt" = "psBayes",
                               "WtRscl" = "psBayes_rscl",
                               "Unwt_sp" = "Bayes_sp",
                               "Wt_sp" = "psBayes_sp",
                               "WtRscl_sp" = "psBayes_rscl_sp",
                               "WtBB" = "psBayes_bb",
                               "WtRsclBB" = "psBayes_rscl_bb")

jittered_locs <- st_jitter(svy_dat, factor = .01)
jittered_locs <- jittered_locs |>
  filter(st_within(jittered_locs, st_as_sf(poly_adm0), sparse = F))

gin_map <- ggplot(data = st_as_sf(poly_adm2)) +
  geom_sf(lwd = .08, fill = NA) + 
  geom_sf(data = st_as_sf(poly_adm1), fill = NA, lwd = .66) + 
  geom_sf(data = jittered_locs %>%
            mutate(urban = as.factor(ifelse(urban == "U", "Urban", "Rural"))),
          aes(color = urban),
          shape = 3, alpha = 1, size = .75) +
  scale_color_manual(values = c("mediumblue", "tomato"), name = NULL) + 
  guides(colour = guide_legend(override.aes = list(size = 4, stroke = 2))) +
  theme_bw() + guides(fill="none") +
  theme(legend.position="bottom",
        legend.text=element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")
ggsave(paste0("paper/figures/gin_map.pdf"),
       width = 8, height = 7)

#### 5.2 MAPS OF ESTIMATES ####################################################
adm2_maps <- st_as_sf(poly_adm2) %>%
  dplyr::select(NAME_2) %>%
  left_join(all_est, by = c("NAME_2" = "GADM"))
all_ests <- ggplot(adm2_maps, aes(fill = median)) + geom_sf(lwd = 0) +
  facet_wrap(~method, ncol = 2) + scale_fill_viridis_c(direction = -1, name = "MCV")  +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) +
  xlab("") + ylab("")
ggsave(paste0("paper/figures/gin_mcv_ests.pdf"),
       width = 10, height = 12)
length_lims <- range(all_est$upper - all_est$lower)
all_lengths <- ggplot(adm2_maps, aes(fill = upper - lower)) +
  geom_sf(lwd = 0) +
  scale_fill_viridis_c(direction = -1, option = "magma",
                       name = "90% CI length") +
  facet_wrap(~method, ncol = 2) + theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())  +
  xlab("") + ylab("")
ggsave(paste0("paper/figures/gin_mcv_int_lengths.pdf"),
       width = 10, height = 12)

library(patchwork)
all_ests + all_lengths
ggsave(paste0("paper/figures/gin_mcv_maps.pdf"),
       width = 20, height = 12)
ht_comp <- all_est %>%
  filter(method != "HAJEK") %>%
  left_join(all_est %>% 
              filter(method == "HAJEK") %>% 
              dplyr::select(domain, est, var) %>% 
              rename(HAJEK = est,
                     HAJEK_var = var),
            by = "domain")

gg <- ggplot(ht_comp, aes(x = HAJEK, y = est, color = method)) +
  geom_point() + geom_abline(slope = 1) + facet_wrap(~method) 
ggsave(plot = gg, height = 6, width = 7,
       filename = paste0("paper/figures/all_adm1_mcv_est_scatter.pdf"))

gg <- ggplot(ht_comp, aes(x = sqrt(HAJEK_var), y = sqrt(var), color = method)) +
  geom_point() + geom_abline(slope = 1) + facet_wrap(~method) 
ggsave(plot = gg, height = 6, width = 7,
       filename = paste0("paper/figures/all_adm1_mcv_se_scatter.pdf"))





#### 6.2 TABLES OF ESTIMATES ####################################################
out_table <- all_est |>
  rename(Prefecture = GADM) |>
  filter(complete.cases(median)) |>
  filter(Prefecture != "Kouroussa") |>
  mutate(est = paste0(round(median, 2), " (", round(lower, 2),", ", round(upper, 2), ")")) |>
  dplyr::select(est, median, method, Prefecture) |>
  filter(method %in% c("HAJEK", "GREG",
                       "Bayes", "psBayes", "psBayes_rscl",
                       "psBayes_bb", "psBayes_rscl_bb")[1:5]) |>
  mutate(method = factor(method, 
                         levels = c("HAJEK", "GREG",
                                    "Bayes", "psBayes", "psBayes_rscl"))) 
levels(out_table$method) <- list("Hájek"  = "HAJEK", 
                                "GREG" = "GREG",
                                "Unwt" = "Bayes",
                                "Wt" = "psBayes",
                                "WtRscl" = "psBayes_rscl") 
out_table <- out_table |>
  pivot_wider(values_from = c("median", "est"),
              names_from = method,
              names_glue = "{method}_{.value}") |>
  arrange(Hájek_median) |>
  dplyr::select(Prefecture, contains("est")) 
names(out_table)[-1] <- stringr::str_replace(names(out_table)[-1], "_est", " (90% CI)")

out_table %>% 
  knitr::kable( format = "latex", booktabs = T, linesep = "") %>%
  writeLines(paste0("paper/figures/gin_mcv_est_table.tex"))
