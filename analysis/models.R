
#### ESTIMATION FUNCTIONS ####
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
                     area_effects = T,
                     CI = .90) {
  
  domain_name <- as.character(domain)[2]
  sample_des$variables$domain_internal <- sample_des$variables[[domain_name]]
  newdata$domain_internal <- newdata[[domain_name]]
  if (area_effects) {
    formula <- update(formula, .~. + domain_internal)
  }
  working_fit <- svyglm(formula,
                        family = family, sample_des)
  pop_unit_ests <- as.vector(predict(working_fit, newdata,
                                     type = "response")) 
  area_ests <-
    aggregate(pop_unit_ests, 
              list(domain_internal = as.vector(newdata$domain_internal)),
              mean)
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
