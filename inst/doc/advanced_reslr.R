## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"#,
  #fig.path = "advanced_vig/"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ---- eval = TRUE, results= 'hide', message=FALSE-----------------------------
#install.packages("reslr")
#devtools::install_github("maeveupton/reslr",force = TRUE,INSTALL_opts = '--no-lock')
library(reslr)

## ----exampledataset2_adv,eval = TRUE------------------------------------------
# For 1 site
CedarIslandNC <- NAACproxydata[NAACproxydata$Site == "Cedar Island",]

## ---- loadigp_adv,eval = TRUE, message=FALSE,results='hide'-------------------
CedarIslandNC_input_detrend <- reslr_load(
  data = CedarIslandNC,
  include_tide_gauge = FALSE,
  include_linear_rate = TRUE,
  TG_minimum_dist_proxy = FALSE,
  list_preferred_TGs = NULL,
  all_TG_1deg = FALSE,
  prediction_grid_res = 50,
  sediment_average_TG = 10,
  detrend_data = TRUE,
  core_col_year = 2010
)

## ----dataigp_adv,eval=FALSE---------------------------------------------------
#  data <- CedarIslandNC_input_detrend$data

## ----datagridigp_adv, eval = FALSE--------------------------------------------
#  data_grid <- CedarIslandNC_input_detrend$data_grid

## ----printigp_adv, eval=TRUE--------------------------------------------------
print(CedarIslandNC_input_detrend)

## ---- plotigpdata_adv,fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(
  x = CedarIslandNC_input_detrend,
  title = "Plot of the raw detrended data",
  xlab = "Year (CE)",
  ylab = "Sea Level (m)",
  plot_proxy_records = TRUE,
  plot_tide_gauges = FALSE
)

## ----runigp_adv,echo = TRUE, eval = FALSE, message=FALSE,results='hide'-------
#  res_eiv_igp_t_detrend <- reslr_mcmc(
#    input_data = CedarIslandNC_input_detrend,
#    model_type = "eiv_igp_t",
#    CI = 0.95
#  )

## ----printigpout_adv, eval=FALSE----------------------------------------------
#  print(res_eiv_igp_t_detrend)

## ----summaryigp_adv, eval = FALSE---------------------------------------------
#  summary(res_eiv_igp_t_detrend)

## ----runigpmore_adv, echo = TRUE, eval = FALSE--------------------------------
#  res_eiv_igp_t_detrend <- reslr_mcmc(
#    input_data = CedarIslandNC_input_detrend,
#    model_type = "eiv_igp_t",
#    # Update these values
#    n_iterations = 6000, # Number of iterations
#    n_burnin = 1000, # Number of iterations to discard at the beginning
#    n_thin = 4, # Reduces number of output samples to save memory and computation time
#    n_chains = 3 # Number of Markov chains
#  )

## ----plotigpres_adv,echo = TRUE, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_eiv_igp_t_detrend,
#    plot_type = "model_fit_plot",
#    xlab = "Year (CE)",
#    ylab = "Sea Level (m)",
#    plot_proxy_records = TRUE,
#    plot_tide_gauges = FALSE
#  )

## ----plotigpresload_adv, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotigpres_adv-1.png"
knitr::include_graphics(url)

## ----plotigpresrate_adv,echo = TRUE, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_eiv_igp_t_detrend,
#    plot_type = "rate_plot",
#    xlab = "Year (CE)",
#    y_rate_lab = "Rate of Change (mm per year)"
#  )

## ----plotigpresrateload_adv, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotigpresrate_adv-1.png"
knitr::include_graphics(url)

## ----igpdataout_adv, eval = FALSE---------------------------------------------
#  output_dataframes <- res_eiv_igp_t_detrend$output_dataframes

## ---- loadBP_adv,eval = FALSE, message=FALSE,results='hide'-------------------
#  CedarIslandNC_input_age_BP <- reslr_load(
#    data = data_age_bp,
#    input_age_type = "BP"
#  )

## ---- eval = FALSE------------------------------------------------------------
#  # For 2 site
#  multi_site <- NAACproxydata[NAACproxydata$Site %in% c("Cedar Island", "Nassau"),]

## ---- eval = FALSE, message=FALSE,results='hide'------------------------------
#  multi_site <- reslr_load(
#    data = multi_site,
#    include_tide_gauge = TRUE,
#    include_linear_rate = TRUE,
#    TG_minimum_dist_proxy = FALSE,
#    # There is no limit to the number of tide gauges provided in the list
#    list_preferred_TGs = c(
#      "ARGENTIA", "MAYPORT",
#      "JACKSONVILLE", "LAKE WORTH PIER",
#      "MAYPORT (BAR PILOTS DOCK), FLORIDA"
#    ),
#    all_TG_1deg = FALSE,
#    prediction_grid_res = 50,
#    sediment_average_TG = 10
#  )

## ----fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE-----------
#  plot(
#    x = multi_site,
#    title = "Plot of the raw data",
#    xlab = "Year (CE)",
#    ylab = "Relative Sea Level (m)",
#    plot_tide_gauges = TRUE,
#    plot_proxy_records = TRUE,
#    plot_caption = TRUE
#  )

## ---- eval = FALSE, message=FALSE,results='hide'------------------------------
#  multi_site <- reslr_load(
#    data = multi_site,
#    include_tide_gauge = TRUE,
#    include_linear_rate = TRUE,
#    TG_minimum_dist_proxy = FALSE,
#    list_preferred_TGs = NULL,
#    all_TG_1deg = TRUE,
#    prediction_grid_res = 50
#  )

## ----fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE-----------
#  plot(
#    x = multi_site,
#    title = "Plot of the raw data",
#    xlab = "Year (CE)",
#    ylab = "Relative Sea Level (m)",
#    plot_tide_gauges = TRUE,
#    plot_proxy_records = TRUE,
#    plot_caption = TRUE
#  )

## ---- eval = FALSE------------------------------------------------------------
#  # Example
#  final_plots <- plot(x = reslr_mcmc(CedarIslandNC, model_type = "ni_spline_t"))
#  final_plots$plot_result
#  # Adding new title to the total model fit plot
#  final_plots$plot_result + ggplot2::ggtitle("New Title Added as Example")
#  final_plots$plot_result + ggplot2::xlab("New x axis label Added as Example")
#  final_plots$plot_result + ggplot2::ylab("New y axis label Added as Example")

## ---- eval=FALSE--------------------------------------------------------------
#  data <- CedarIslandNC_input_detrend$data

## ---- eval=FALSE--------------------------------------------------------------
#  data <- res_eiv_igp_t_detrend$output_dataframes

## ---- eval = FALSE------------------------------------------------------------
#  # Example
#  CedarIslandNC_input <- reslr_load(
#    data = CedarIslandNC)
#  res_eiv_slr_t <-
#    reslr_mcmc(CedarIslandNC_input,
#               model_type = "eiv_slr_t")
#  # Accessing the slope of the EIV simple linear regression
#  beta <- res_eiv_slr_t$noisy_model_run_output$BUGSoutput$sims.list$beta

## ---- eval = FALSE------------------------------------------------------------
#  res_ni_sp_t <-
#    reslr_mcmc(CedarIslandNC_input,
#               model_type = "ni_spline_t",
#               spline_nseg = NULL)

## ---- eval = FALSE------------------------------------------------------------
#  res_ni_sp_t <-
#    reslr_mcmc(CedarIslandNC_input,
#               model_type = "ni_spline_st",
#               spline_nseg = NULL)

## ---- eval = FALSE------------------------------------------------------------
#  res_ni_sp_t <-
#    reslr_mcmc(CedarIslandNC_input,
#               model_type = "ni_gam_decomp",
#               spline_nseg_t = 20,
#               spline_nseg_st = 6)

## ----cv_adv,eval = TRUE,message=FALSE,results='hide'--------------------------
data1site_example <- NAACproxydata[NAACproxydata$Site == "Cedar Island",]
# Cross Validation test
cv <- cross_val_check(data = data1site_example,
                      model_type ="ni_spline_t",
                      n_iterations = 1000,
                      n_burnin = 100,
                      n_thin = 5,
                      n_chains = 2,
                      spline_nseg = NULL,# User the package to calculate the number of knots
                      # n_fold allows the user to alter the cross validation, i.e. 3, 5, 10 fold
                      n_fold = 3,
                      #To reproducible results,seed stores the output of the random selection
                      seed = NULL,
                      CI = 0.95)# Size of the credible intervals and prediction intervals


## ----cvpredplot_adv,eval = TRUE-----------------------------------------------
cv$true_pred_plot

## ----cvdf_adv,eval = FALSE----------------------------------------------------
#  CV_model_df <- cv$CV_model_df

## ----coverage_adv,eval = TRUE-------------------------------------------------
# Overall coverage
total_empirical_coverage <- cv$total_coverage
total_empirical_coverage
# Coverage by site
coverage_by_site <- cv$coverage_by_site
coverage_by_site
# Size of the prediction intervals
prediction_interval_size <- cv$prediction_interval_size
prediction_interval_size

## ----rmse_adv,eval = TRUE-----------------------------------------------------
# Overall
ME_MAE_RSME_overall <- cv$ME_MAE_RSME_overall
ME_MAE_RSME_overall
# By fold and site
ME_MAE_RSME_fold_site <- cv$ME_MAE_RSME_fold_site
ME_MAE_RSME_fold_site
# By site
ME_MAE_RSME_site <- cv$ME_MAE_RSME_site
ME_MAE_RSME_site
# By fold
ME_MAE_RSME_fold <- cv$ME_MAE_RSME_fold
ME_MAE_RSME_fold

