library(data.table)

f <- function(y) {
  last <- !duplicated(y$t, fromLast = TRUE)
  last_ind <- which(last == TRUE)
  return(seq(0, y$t[last_ind]))
}

trial_period_func <- function(x) {
  # Dummy variables used in data.table calls declared to prevent package check NOTES:
  t <- ID <- trial_period <- NULL
  
  x_new <- x[rep(seq_len(.N), t + 1), list(ID, t)]
  x_new[, trial_period := f(.BY), by = list(ID, t)]
  return(x_new[, trial_period])
}

quiet_print <- function(quiet, x, ...) {
  if (isFALSE(quiet)) {
    print(x, ...)
  }
}

#' Conditional Messages
#'
#' @param quiet (`logical`) Messages printed if `FALSE`
#' @param x Object to print.
#' @param ... Passed to `message` function.
#'
#' @noRd
quiet_msg <- function(quiet, x, ...) {
  if (isFALSE(quiet)) {
    message(x, ...)
  }
}

#' Print a line
#'
#' @param quiet Print if `TRUE `
#' @noRd
quiet_line <- function(quiet) {
  quiet_msg(quiet, paste0(strrep("-", 0.75 * getOption("width")), "\n"))
}

#' Print with timing statement
#'
#' @param quiet Print if `TRUE `,
#' @param x Message to print
#' @param proc_time Result of `system.time()`. Elapsed time will be extracted,
#' formatted for printing and `paste0()`ed to `x`.
#' @noRd
quiet_msg_time <- function(quiet, msg, proc_time) {
  time <- proc_time["elapsed"]
  time <- if (time < 10) sprintf("%0.1f s", time) else sprintf("%.5g s", time)
  quiet_msg(quiet, paste0(msg, time))
}

#' Weight model fitting on existing expanded data
#'
#' @param quiet Print if `TRUE '
#' @param data
#' @param expanded_data Expanded data, output of data_preparation
#' @param switch_d_cov model formula for switching weight model at the denominator
#' @param switch_n_cov 
#' @param cense 
#' @param pool_cense
#' @param cense_d_cov 
#' @param cense_n_cov  
#' @param eligible_wts_0 
#' @param eligible_wts_1 
#' @param remodel TRUE by default, indicator of whether weight models should be fitted/refitted
#' @param boot_idx 2-column matrix of patient IDs and their bootstrap resampling weight. I.e. if the patient is sammpled 2 times in the 
#' bootstrap sample, weight is 2. If they weren't sampled, weight is 0.
#' @param include_regime_length 
#' @param weight_model_d0_data 
#' @param weight_model_n0_data 
#' @param weight_model_d1_data 
#' 
#' @noRd
weight_func_bootstrap <- function(data,
                                  expanded_data,
                                  treatment = 'A',
                                  switch_n_cov = ~ 1,
                                  switch_d_cov = NA,
                                  cense = NA,
                                  pool_cense = 0,
                                  cense_d_cov = NA,
                                  cense_n_cov = ~ 1,
                                  eligible_wts_0 = NA,
                                  eligible_wts_1 = NA,
                                  include_regime_length = 0,
                                  weight_model_d0 = switch_d0,
                                  weight_model_n0 = switch_n0,
                                  weight_model_d1 = switch_d1,
                                  weight_model_n1 = switch_n1,
                                  cense_model_d0 = cense_d0,
                                  cense_model_n0 = cense_n0,
                                  cense_model_d1 = cense_d1,
                                  cense_model_n1 = cense_n1,
                                  remodel = TRUE,
                                  new_coef_sw_d0 = NA,
                                  new_coef_sw_n0 = NA,
                                  new_coef_sw_d1 = NA,
                                  new_coef_sw_n1 = NA,
                                  new_coef_c_d0 = NA,
                                  new_coef_c_n0 = NA,
                                  new_coef_c_d1 = NA,
                                  new_coef_c_n1 = NA,
                                  boot_idx = NA,
                                  save_weight_models = FALSE,
                                  quiet = FALSE) {
  
  
  if (all(!is.na(boot_idx), na.rm = TRUE)){ #Modify weight model design matrix to only include sampled IDs if boot_idx is not NA
    weight_model_d0_data <- as.data.table(weight_model_d0$data[id %in% boot_idx] %>% 
      rowwise() %>% 
      dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
    weight_model_n0_data <- as.data.table(weight_model_n0$data[id %in% boot_idx] %>% 
      rowwise() %>% 
      dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
    weight_model_d1_data <- as.data.table(weight_model_d1$data[id %in% boot_idx]%>% 
      rowwise() %>% 
      dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
    weight_model_n1_data <- as.data.table(weight_model_n1$data[id %in% boot_idx]%>% 
      rowwise() %>% 
      dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
    if(!is.na(cense)){
      cense_model_d0_data <- as.data.table(cense_model_d0$data[id %in% boot_idx]%>% 
        rowwise() %>% 
        dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
      cense_model_n0_data <- as.data.table(cense_model_n0$data[id %in% boot_idx]%>% 
        rowwise() %>% 
        dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
      cense_model_d1_data <- as.data.table(cense_model_d1$data[id %in% boot_idx]%>% 
        rowwise() %>% 
        dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
      cense_model_n1_data <- as.data.table(cense_model_n1$data[id %in% boot_idx]%>% 
        rowwise() %>% 
        dplyr::mutate(weight_boot = length(boot_idx[boot_idx == id])))
    }
  } else {
    weight_model_d0_data <- as.data.table(weight_model_d0$data %>% 
      dplyr::mutate(weight_boot = 1.0))
    weight_model_n0_data <- as.data.table(weight_model_n0$data %>% 
      dplyr::mutate(weight_boot = 1.0))
    weight_model_d1_data <- as.data.table(weight_model_d1$data %>% 
      dplyr::mutate(weight_boot = 1.0))
    weight_model_n1_data <- as.data.table(weight_model_n1$data %>% 
      dplyr::mutate(weight_boot = 1.0))
    if (!is.na(cense)){
      cense_model_d0_data <- as.data.table(cense_model_d0$data %>% 
        dplyr::mutate(weight_boot = 1.0))
      cense_model_n0_data <- as.data.table(cense_model_n0$data %>% 
        dplyr::mutate(weight_boot = 1.0))
      cense_model_d1_data <- as.data.table(cense_model_d1$data %>% 
        dplyr::mutate(weight_boot = 1.0))
      cense_model_n1_data <- as.data.table(cense_model_n1$data %>% 
        dplyr::mutate(weight_boot = 1.0))
    }
  }
  
  eligible0 <- eligible1 <- ID <- t <- eligible0.y <- eligible1.y <- am_1 <-
    treatment <- wt <- wtC <- p0_n <- p0_d <- p1_n <- p1_d <- pC_n0 <- pC_d0 <-
    pC_n1 <- pC_d1 <- pC_n <- pC_d <- NULL
  
  data <- data %>% 
    mutate(eligible0 = ifelse(Ap == 0, 1, 0),
           eligible1 = ifelse(Ap == 1, 1 , 0))
  
  switch_d_cov <- update.formula(switch_d_cov, treatment ~ .)
  switch_n_cov <- update.formula(switch_n_cov, treatment ~ .)
  
  if (include_regime_length == 1) {
    switch_d_cov <- update.formula(switch_d_cov, ~ . + time_on_regime + I(time_on_regime^2))
    switch_n_cov <- update.formula(switch_n_cov, ~ . + time_on_regime + I(time_on_regime^2))
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Switching weights --------------------
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # Fit the models for the weights in the four scenarios
  weight_models <- list()
  # ------------------- eligible0 == 1 --------------------
  # --------------- denominator ------------------
  quiet_msg(quiet, "P(treatment=1 | treatment=0) for denominator")
  if (remodel ==TRUE){ #fit/refit weight models 
    model1 <- parglm::parglm(switch_d_cov,
                             data = weight_model_d0_data,
                             weights = weight_boot,
                             family = binomial(link = "logit"),
                             control = parglm::parglm.control(nthreads = 4, method = "FAST"))
    quiet_print(quiet, summary(model1))
    switch_d0 <- data.table(
      p0_d = model1$fitted.values,
      eligible0 = unlist(model1$data$eligible0),
      ID = model1$data[, id],
      t = model1$data[, period]
    )
    weight_models$switch_d0 <- broom::tidy(model1)
    weight_models$switch_d0_statistics <- broom::glance(model1)
    if (save_weight_models) {
      saveRDS(model1, file = file.path(save_dir, "weight_model_switch_d0.rds"))
    }
    rm(model1)
  } else { #if not remodelling, predict probability on new data/bootstrap sample
    weight_model_d0$coefficients <- new_coef_sw_d0
    switch_d0 <- data.table(
      p0_d = predict.glm(weight_model_d0, weight_model_d0_data, type = 'response' ),
      eligible0 = unlist(weight_model_d0_data$eligible0),
      ID = weight_model_d0_data[, id],
      t = weight_model_d0_data[, period])
  }
  
  # -------------- numerator --------------------
  
  quiet_msg(quiet, "P(treatment=1 | treatment=0) for numerator")
  if (remodel ==TRUE){
    model2 <- parglm::parglm(switch_n_cov,
                             data = weight_model_n0_data,
                             weights = weight_boot,
                             family = binomial(link = "logit"),
                             control = parglm::parglm.control(nthreads = 4, method = "FAST"))
    quiet_print(quiet, summary(model2))
    switch_n0 <- data.table(
      p0_n = model2$fitted.values,
      eligible0 = unlist(model2$data$eligible0),
      ID = model2$data[, id],
      t = model2$data[, period]
    )
    weight_models$switch_n0 <- broom::tidy(model2)
    weight_models$switch_n0_statistics <- broom::glance(model2)
    if (save_weight_models) {
      saveRDS(model2, file = file.path(save_dir, "weight_model_switch_n0.rds"))
    }
    rm(model2)
  } else {
    weight_model_n0$coefficients <- new_coef_sw_n0
    switch_n0 <- data.table(
      p0_n = predict.glm(weight_model_n0, weight_model_n0_data, type = 'response' ),
      eligible0 = unlist(weight_model_n0_data$eligible0),
      ID = weight_model_n0_data[, id],
      t = weight_model_n0_data[, period])
  }
  
  # ------------------- eligible1 == 1 --------------------
  # --------------- denominator ------------------
  quiet_msg(quiet, "P(treatment=1 | treatment=1) for denominator")
  if (remodel ==TRUE){
    model3 <- parglm::parglm(switch_d_cov,
                             data = weight_model_d1_data,
                             weights = weight_boot,
                             family = binomial(link = "logit"),
                             control = parglm::parglm.control(nthreads = 4, method = "FAST"))
    quiet_print(quiet, summary(model3))
    switch_d1 <- data.table(
      p1_d = model3$fitted.values,
      eligible1 = unlist(model3$data$eligible1),
      ID = model3$data[, id],
      t = model3$data[, period]
    )
    weight_models$switch_d1 <- broom::tidy(model3)
    weight_models$switch_d1_statistics <- broom::glance(model3)
    if (save_weight_models) {
      saveRDS(model3, file = file.path(save_dir, "weight_model_switch_d1.rds"))
    }
    rm(model3)
  } else {
    weight_model_d1$coefficients <- new_coef_sw_d1
    switch_d1 <- data.table(
      p1_d = predict.glm(weight_model_d1, weight_model_d1_data, type = 'response' ),
      eligible1 = unlist(weight_model_d1_data$eligible1),
      ID = weight_model_d1_data[, id],
      t = weight_model_d1_data[, period])
  }
  
  # -------------------- numerator ---------------------------
  quiet_msg(quiet, "P(treatment=1 | treatment=1) for numerator")
  if (remodel ==TRUE){
    model4 <- parglm::parglm(switch_n_cov,
                             data = weight_model_n1_data,
                             weights = weight_boot,
                             family = binomial(link = "logit"),
                             control = parglm::parglm.control(nthreads = 4, method = "FAST"))
    quiet_print(quiet, summary(model4))
    switch_n1 <- data.table(
      p1_n = model4$fitted.values,
      eligible1 = unlist(model4$data$eligible1),
      ID = model4$data[, id],
      t = model4$data[, period]
    )
    weight_models$switch_n1 <- broom::tidy(model4)
    weight_models$switch_n1_statistics <- broom::glance(model4)
    if (save_weight_models) {
      saveRDS(model4, file = file.path(save_dir, "weight_model_switch_n1.rds"))
    }
    rm(model4)
  } else {
    weight_model_n1$coefficients <- new_coef_sw_n1
    switch_n1 <- data.table(
      p1_n = predict.glm(weight_model_n1, weight_model_n1_data, type = 'response' ),
      eligible1 = unlist(weight_model_n1_data$eligible1),
      ID = weight_model_n1_data[, id],
      t = weight_model_n1_data[, period])
  }
  
  print('models refitted')
  # -------------- Combine results --------------------
  
  switch_0 <- switch_d0[switch_n0, on = list(
    ID = ID, t = t,
    eligible0 = eligible0
  )]
  switch_1 <- switch_d1[switch_n1, on = list(
    ID = ID, t = t,
    eligible1 = eligible1
  )]
  
  rm(switch_d0, switch_d1, switch_n0, switch_n1)
  
  data <- merge.data.table(data, switch_0, by = c("ID", "t"), all = TRUE)
  data <- merge.data.table(data, switch_1, by = c("ID", "t"), all = TRUE)
  
  rm(switch_1, switch_0)
  data <- data.table(data) 
  data[, eligible0.y := NULL]
  data[, eligible1.y := NULL]
  setnames(data, c("eligible0.x", "eligible1.x"), c("eligible0", "eligible1"))
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Censoring weights --------------------
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cens_models <- list()
  if (!is.na(cense)) {
    cense_d_cov <- update(cense_d_cov, paste("1 -", cense, "~ ."))
    cense_n_cov <- update(cense_n_cov, paste("1 -", cense, "~ ."))
    
    if (pool_cense == 1) {
      # -------------------- denominator -------------------------
      quiet_msg(quiet, "Model for P(cense = 0 |  X ) for denominator")
      # -----------------------------------------------------------
      model1.cense <- parglm::parglm(cense_d_cov,
                                     data = data,
                                     weights = weight_boot,
                                     family = binomial(link = "logit"),
                                     control = parglm::parglm.control(nthreads = 4, method = "FAST"))
      quiet_print(quiet, summary(model1.cense))
      cense_d0 <- data.table(
        pC_d = model1.cense$fitted.values,
        ID = model1.cense$data[, id],
        t = model1.cense$data[, period]
      )
      
      cens_models$cens_pool_d <- broom::tidy(model1.cense)
      cens_models$cens_pool_d_statistics <- broom::glance(model1.cense)
      if (save_weight_models) {
        saveRDS(model1.cense, file = file.path(save_dir, "cense_model_pool_d.rds"))
      }
      rm(model1.cense)
      
      # --------------------- numerator ---------------------------
      quiet_msg(quiet, "Model for P(cense = 0 |  X ) for numerator")
      # ---------------------------------------------------------
      model2.cense <- parglm::parglm(cense_n_cov,
                                     data = data,
                                     weights = weight_boot,
                                     family = binomial(link = "logit"),
                                     control = parglm::parglm.control(nthreads = 4, method = "FAST"))
      quiet_print(quiet, summary(model2.cense))
      cense_n0 <- data.table(
        pC_n = model2.cense$fitted.values,
        ID = model2.cense$data[, id],
        t = model2.cense$data[, period]
      )
      
      cens_models$cens_pool_n <- broom::tidy(model2.cense)
      cens_models$cens_pool_n_statistics <- broom::glance(model2.cense)
      if (save_weight_models) {
        saveRDS(model2.cense, file = file.path(save_dir, "cense_model_pool_n.rds"))
      }
      rm(model2.cense)
      data <- merge.data.table(data, cense_d0, by = c("ID", "t"), all = TRUE)
      data <- merge.data.table(data, cense_n0, by = c("ID", "t"), all = TRUE)
      
      rm(cense_d0, cense_n0)
    } else {
      # when pool_cense != 1
      
      # ---------------------- denominator -----------------------
      quiet_msg(quiet, "Model for P(cense = 0 |  X, Am1=0) for denominator")
      # ---------------------- eligible0 ---------------------------
      if (remodel == TRUE){
        model1.cense <- parglm::parglm(cense_d_cov,
                                       data = cense_model_d0_data,
                                       weights = weight_boot,
                                       family = binomial(link = "logit"),
                                       control = parglm::parglm.control(nthreads = 4, method = "FAST"))
        quiet_print(quiet, summary(model1.cense))
        cense_d0 <- data.table(
          pC_d0 = model1.cense$fitted.values,
          ID = model1.cense$data[, id],
          t = model1.cense$data[, period]
        )
        cens_models$cens_d0 <- broom::tidy(model1.cense)
        cens_models$cens_d0_statistics <- broom::glance(model1.cense)
        if (save_weight_models) {
          saveRDS(model1.cense, file = file.path(save_dir, "cense_model_d0.rds"))
        }
        rm(model1.cense)
      } else {
        cense_model_d0$coefficients <- new_coef_c_d0
        cense_d0 <- data.table(
          pC_d0 = predict.glm(cense_model_d0, cense_model_d0_data, type = 'response'),
          ID = cense_model_d0_data[, id],
          t = cense_model_d0_data[, period]
        )
      }
      
      
      # -------------------------- numerator ----------------------
      quiet_msg(quiet, "Model for P(cense = 0 |  X, Am1=0) for numerator")
      #--------------------------- eligible0 -----------------------
      if (remodel == TRUE){
        model2.cense <- parglm::parglm(cense_n_cov,
                                       data = cense_model_n0_data,
                                       weights = weight_boot,
                                       family = binomial(link = "logit"),
                                       control = parglm::parglm.control(nthreads = 4, method = "FAST"))
        quiet_print(quiet, summary(model2.cense))
        cense_n0 <- data.table(
          pC_n0 = model2.cense$fitted.values,
          ID = model2.cense$data[, id],
          t = model2.cense$data[, period]
        )
        cens_models$cens_n0 <- broom::tidy(model2.cense)
        cens_models$cens_n0_statistics <- broom::glance(model2.cense)
        if (save_weight_models) {
          saveRDS(model2.cense, file = file.path(save_dir, "cense_model_n0.rds"))
        }
        rm(model2.cense)
      } else {
        cense_model_n0$coefficients <- new_coef_c_n0
        cense_n0 <- data.table(
          pC_n0 = predict.glm(cense_model_n0, cense_model_n0_data, type = 'response'),
          ID = cense_model_n0_data[, id],
          t = cense_model_n0_data[, period]
        )
      }
      
      # ------------------------- denominator ---------------------
      quiet_msg(quiet, "Model for P(cense = 0 |  X, Am1=1) for denominator")
      # ------------------------ eligible1 -------------------------
      if (remodel == TRUE){
        model3.cense <- parglm::parglm(cense_d_cov,
                                       data = cense_model_d1_data,
                                       weights = weight_boot,
                                       family = binomial(link = "logit"),
                                       control = parglm::parglm.control(nthreads = 4, method = "FAST"))
        quiet_print(quiet, summary(model3.cense))
        cense_d1 <- data.table(
          pC_d1 = model3.cense$fitted.values,
          ID = model3.cense$data[, id],
          t = model3.cense$data[, period]
        )
        cens_models$cens_d1 <- broom::tidy(model3.cense)
        cens_models$cens_d1_statistics <- broom::glance(model3.cense)
        if (save_weight_models) {
          saveRDS(model3.cense, file = file.path(save_dir, "cense_model_d1.rds"))
        }
        rm(model3.cense)
      } else {
        cense_model_d1$coefficients <- new_coef_c_d1
        cense_d1 <- data.table(
          pC_d1 = predict.glm(cense_model_d1, cense_model_d1_data, type = 'response'),
          ID = cense_model_d1_data[, id],
          t = cense_model_d1_data[, period]
        )
      }
      # ------------------------ numerator -------------------------
      quiet_msg(quiet, "Model for P(cense = 0 |  X, Am1=1) for numerator")
      # ------------------------- eligible1 -----------------------
      if (remodel == TRUE){
        model4.cense <- parglm::parglm(cense_n_cov,
                                       data = cense_model_n1_data,
                                       weights = weight_boot,
                                       family = binomial(link = "logit"),
                                       control = parglm::parglm.control(nthreads = 4, method = "FAST"))
        quiet_print(quiet, summary(model4.cense))
        cense_n1 <- data.table(
          pC_n1 = model4.cense$fitted.values,
          ID = model4.cense$data[, id],
          t = model4.cense$data[, period]
        )
        cens_models$cens_n1 <- broom::tidy(model4.cense)
        cens_models$cens_n1_statistics <- broom::glance(model4.cense)
        if (save_weight_models) {
          saveRDS(model4.cense, file = file.path(save_dir, "cense_model_n1.rds"))
        }
        rm(model4.cense)
      } else {
        cense_model_n1$coefficients <- new_coef_c_n1
        cense_n1 <- data.table(
          pC_n1 = predict.glm(cense_model_n1, cense_model_n1_data, type = 'response'),
          ID = cense_model_n1_data[, id],
          t = cense_model_n1_data[, period]
        )
      }
      
      # combine ------------------------------
      cense_0 <- cense_d0[cense_n0, on = list(ID = ID, t = t)]
      cense_1 <- cense_d1[cense_n1, on = list(ID = ID, t = t)]
      rm(cense_n1, cense_d1, cense_n0, cense_d0)
      data <- merge.data.table(data, cense_0, by = c("ID", "t"), all = TRUE)
      data <- merge.data.table(data, cense_1, by = c("ID", "t"), all = TRUE)
      rm(cense_0, cense_1)
    }
  }
  
  quiet_msg(quiet, "Calculating weights")
  # wt and wtC calculation
  if (any(!is.na(eligible_wts_0))) {
    data[
      ( Ap == 0 & eligible == 1 & A == 0 & !is.na(p0_n) & !is.na(p0_d)),
      wt := (1.0 - p0_n) / (1.0 - p0_d)
    ]
    data[
      (Ap == 0 & eligible == 1 & A == 1 & !is.na(p0_n) & !is.na(p0_d)),
      wt := p0_n / p0_d
    ]
    data[(Ap == 0 & eligible == 0), wt := 1.0]
  } else {
    data[
      (Ap == 0 & A == 0 & !is.na(p0_n) & !is.na(p0_d)),
      wt := (1.0 - p0_n) / (1.0 - p0_d)
    ]
    data[
      (Ap == 0 & A == 1 & !is.na(p0_n) & !is.na(p0_d)),
      wt := p0_n / p0_d
    ]
  }
  if (any(!is.na(eligible_wts_1))) {
    data[
      (Ap == 1 & eligible == 1 & A == 0 & !is.na(p1_n) & !is.na(p1_d)),
      wt := (1.0 - p1_n) / (1.0 - p1_d)
    ]
    data[
      (Ap == 1 & eligible == 1 & A == 1 & !is.na(p1_n) & !is.na(p1_d)),
      wt := p1_n / p1_d
    ]
    data[(Ap == 1 & eligible == 0), wt := 1.0]
  } else {
    data[
      (Ap == 1 & A == 0 & !is.na(p1_n) & !is.na(p1_d)),
      wt := (1.0 - p1_n) / (1.0 - p1_d)
    ]
    data[
      (Ap == 1 & A == 1 & !is.na(p1_n) & !is.na(p1_d)),
      wt := p1_n / p1_d
    ]
  }
  
  if (is.na(cense)) {
    data[, wtC := 1.0]
  } else {
    if (pool_cense == 0) {
      data[Ap == 0, `:=`(pC_n = pC_n0, pC_d = pC_d0)]
      data[Ap == 1, `:=`(pC_n = pC_n1, pC_d = pC_d1)]
    }
    data[is.na(pC_d), pC_d := 1]
    data[is.na(pC_n), pC_n := 1]
    data[, wtC := pC_n / pC_d]
  }
  data[, wt := wt * wtC]
  data[, first := !duplicated(data[,ID])]
  data <- data[!is.na(wt)]
  temp_data <- data.table(
    id = data[, ID],
    period = data[, t]
  )
  temp_data[, wtprod := 1.0, by = id]
  #[, elgcount := 0.0, by = id][, expand := 0.0, by = id]
  #temp_data[, treat := 0.0, by = id][, dosesum := 0.0, by = id]
  
  data[first == TRUE, weight0 := 1.0]
  data[, weight0 := cumprod(wt), by = ID]
  temp_data[, wtprod := data[, weight0]]
  #temp_data[, treat := data[, A]]
  #temp_data[, dosesum := data[, CAp]]
  #temp_data[, elgcount := data[, eligible]]
  #temp_data[data[, eligible] == 1, init := data[eligible == 1, A]]
  #temp_data[, init_shift := shift(data[, A])]
  #temp_data[data[, eligible] == 0, init := init_shift, by = id]
  #temp_data[, init_shift := NULL]
  
  expand_index <- rep(seq_len(nrow(data)), data[, t] + 1)
  
  quiet_msg(quiet, "Adding new weights to expanded data")
  
  ### new_data only contains ID, trial_period, followup_time adn the new IP weights
  new_data <- data.table(id = data[expand_index, ID])
  new_data[, period_new := data[expand_index, t]]
  #new_data[, cumA_new := data[expand_index, CAp]]
  #new_data[, treatment_new := data[expand_index, A]]
  
  quiet_msg(quiet, "Placer 1")
  #new_data[, outcome_new := data[expand_index, Y]]
  new_data[, weight0 := data[expand_index, weight0]]
  new_data[, trial_period := trial_period_func(data)]
  new_data[, index := seq_len(.N)]
  
  quiet_msg(quiet, "Placer 2")
  new_data <- new_data[temp_data, on = list(id = id, trial_period = period)]
  setorder(new_data, index)
  new_data[, followup_time := period_new - trial_period]
  new_data[, weight := (weight0 / wtprod)]
  
  quiet_msg(quiet, "Placer 3")
  #### New data is merged with existing expanded data to add the new weights 
  output_data <- new_data[expanded_data, on = list( id = id, trial_period = trial_period, 
                                                    followup_time = followup_time)] %>%
    dplyr::select(names(expanded_data))
  
  list(
    data = output_data,
    switch_models = weight_models,
    censor_models = cens_models
  )
}
