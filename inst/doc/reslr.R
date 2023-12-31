## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"#,
  #fig.path = "vig/"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----eval = FALSE,message=FALSE-----------------------------------------------
#  # install.packages("reslr")
#  # library(devtools)
#  # devtools::install()
#  #devtools::install_github("maeveupton/reslr")
#  install_github("maeveupton/reslr")

## ---- readpkg,eval = TRUE, message=FALSE--------------------------------------
library(reslr)

## ---- example1site,eval = TRUE------------------------------------------------
# For 1 site
data_1site <- reslr::NAACproxydata %>% dplyr::filter(Site == "Cedar Island")
# For multiple sites
data_multisite <- reslr::NAACproxydata %>% dplyr::filter(Site %in% c(
  "Snipe Key", "Cheesequake",
  "Placentia", "Leeds Point"
))

## ---- example1dataagain, eval = TRUE------------------------------------------
# For 1 site
CedarIslandNC <- NAACproxydata %>% dplyr::filter(Site == "Cedar Island")

## ---- loadslr, eval = TRUE----------------------------------------------------
CedarIslandNC_input <- reslr_load(
  data = CedarIslandNC,
  include_tide_gauge = FALSE,
  include_linear_rate = FALSE,
  TG_minimum_dist_proxy = FALSE,
  list_preferred_TGs = NULL,
  all_TG_1deg = FALSE,
  prediction_grid_res = 50,
  input_age_type = "CE",
  sediment_average_TG = 10
)

## ---- slrdata,eval=TRUE-------------------------------------------------------
data <- CedarIslandNC_input$data

## ---- datagridslr, eval = TRUE------------------------------------------------
data_grid <- CedarIslandNC_input$data_grid

## ---- printdata, eval=TRUE----------------------------------------------------
print(CedarIslandNC_input)

## ---- plotdata,fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(
  x = CedarIslandNC_input,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_tide_gauges = FALSE,
  plot_proxy_records = TRUE,
  plot_caption = TRUE
)

## ---- runslr,eval = TRUE------------------------------------------------------
res_eiv_slr_t <- reslr_mcmc(
  input_data = CedarIslandNC_input,
  model_type = "eiv_slr_t",
  CI = 0.95
)

## ---- printrunslr,eval=TRUE---------------------------------------------------
print(res_eiv_slr_t)

## ---- summaryslr, eval = TRUE-------------------------------------------------
summary(res_eiv_slr_t)

## ---- runslrmoreit,eval = FALSE-----------------------------------------------
#  res_eiv_slr_t <- reslr_mcmc(
#    input_data = CedarIslandNC_input,
#    model_type = "eiv_slr_t",
#    # Update these values
#    n_iterations = 6000, # Number of iterations
#    n_burnin = 1000, # Number of iterations to discard at the beginning
#    n_thin = 4, # Reduces number of output samples to save memory and computation time
#    n_chains = 3 # Number of Markov chains
#  )

## ---- plotslrres, fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(res_eiv_slr_t,
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)"
)

## ---- dataframeslrres, eval = TRUE--------------------------------------------
output_dataframes <- res_eiv_slr_t$output_dataframes
head(output_dataframes)

## ----exampledata1, eval = TRUE------------------------------------------------
# For 1 site
CedarIslandNC <- reslr::NAACproxydata %>% dplyr::filter(Site == "Cedar Island")

## ----loaddatacp, eval = TRUE--------------------------------------------------
CedarIslandNC_input <- reslr_load(
  data = CedarIslandNC,
  include_tide_gauge = FALSE,
  include_linear_rate = FALSE,
  TG_minimum_dist_proxy = FALSE,
  list_preferred_TGs = NULL,
  all_TG_1deg = FALSE,
  prediction_grid_res = 50,
  sediment_average_TG = 10
)

## ----datacp,eval=TRUE---------------------------------------------------------
data <- CedarIslandNC_input$data
head(data)

## ----datagridcp, eval = TRUE--------------------------------------------------
data_grid <- CedarIslandNC_input$data_grid
head(data_grid)

## ----printcp, eval=TRUE-------------------------------------------------------
print(CedarIslandNC_input)

## ----plotdatacp,fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(
  x = CedarIslandNC_input,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_proxy_records = TRUE,
  plot_tide_gauges = FALSE
)

## ----runcp1, eval = TRUE------------------------------------------------------
res_eiv_cp1_t <- reslr_mcmc(
  input_data = CedarIslandNC_input,
  model_type = "eiv_cp_t",
  n_cp = 1,
  CI = 0.95
)

## ----runcp2, eval = FALSE-----------------------------------------------------
#  res_eiv_cp2_t <- reslr_mcmc(
#    input_data = CedarIslandNC_input,
#    model_type = "eiv_cp_t",
#    n_cp = 2, # Updating the default setting to include an additional change point.
#    CI = 0.95
#  )

## ----printcpout, eval=TRUE----------------------------------------------------
print(res_eiv_cp1_t)

## ----summarycp, eval = TRUE---------------------------------------------------
summary(res_eiv_cp1_t)

## ----runcpmoreit,eval = FALSE-------------------------------------------------
#  res_eiv_cp1_t <- reslr_mcmc(
#    input_data = CedarIslandNC_input,
#    model_type = "eiv_cp_t",
#    # Update these values
#    n_iterations = 6000, # Number of iterations
#    n_burnin = 1000, # Number of iterations to discard at the beginning
#    n_thin = 4, # Reduces number of output samples to save memory and computation time
#    n_chains = 3 # Number of Markov chains
#  )

## ---- plotcpres, fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(res_eiv_cp1_t,
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
)

## ----cpdataframe, eval = TRUE-------------------------------------------------
output_dataframes <- res_eiv_cp1_t$output_dataframes
head(output_dataframes)

## ----exampledataset2 ,eval = TRUE---------------------------------------------
# For 1 site
CedarIslandNC <- reslr::NAACproxydata %>% dplyr::filter(Site == "Cedar Island")

## ---- loadigp ,eval = TRUE----------------------------------------------------
CedarIslandNC_input <- reslr_load(
  data = CedarIslandNC,
  include_tide_gauge = FALSE,
  include_linear_rate = FALSE,
  TG_minimum_dist_proxy = FALSE,
  list_preferred_TGs = NULL,
  all_TG_1deg = FALSE,
  prediction_grid_res = 50,
  sediment_average_TG = 10
)

## ----dataigp,eval=FALSE-------------------------------------------------------
#  data <- CedarIslandNC_input$data

## ----datagridigp, eval = FALSE------------------------------------------------
#  data_grid <- CedarIslandNC_input$data_grid

## ----printigp, eval=TRUE------------------------------------------------------
print(CedarIslandNC_input)

## ---- plotigpdata,fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(
  x = CedarIslandNC_input,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_proxy_records = TRUE,
  plot_tide_gauges = FALSE
)

## ----runigp, eval = FALSE-----------------------------------------------------
#  res_eiv_igp_t <- reslr_mcmc(
#    input_data = CedarIslandNC_input,
#    model_type = "eiv_igp_t",
#    CI = 0.95,
#  
#  )

## ----printigpout, eval=FALSE--------------------------------------------------
#  print(res_eiv_igp_t)

## ----summaryigp, eval = FALSE-------------------------------------------------
#  summary(res_eiv_igp_t)

## ----runigpmore, eval = FALSE-------------------------------------------------
#  res_eiv_igp_t <- reslr_mcmc(
#    input_data = CedarIslandNC_input,
#    model_type = "eiv_igp_t",
#    # Update these values
#    n_iterations = 6000, # Number of iterations
#    n_burnin = 1000, # Number of iterations to discard at the beginning
#    n_thin = 4, # Reduces number of output samples to save memory and computation time
#    n_chains = 3 # Number of Markov chains
#  )

## ----plotigpres, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_eiv_igp_t,
#    plot_type = "model_fit_plot",
#    xlab = "Year (CE)",
#    ylab = "Relative Sea Level (m)",
#    plot_proxy_records = TRUE,
#    plot_tide_gauges = FALSE
#  )

## ----plotigpres_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotigpres-1.png"
knitr::include_graphics(url)

## ----plotigpresrate, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_eiv_igp_t,
#    plot_type = "rate_plot",
#    xlab = "Year (CE)",
#    y_rate_lab = "Rate of Change (mm per year)"
#  )

## ----plotigpresrate_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotigpresrate-1.png"
knitr::include_graphics(url)

## ----igpdataout, eval = FALSE-------------------------------------------------
#  output_dataframes <- res_eiv_igp_t$output_dataframes

## ---- dataexample3,eval = TRUE------------------------------------------------
# For 1 site
CedarIslandNC <- reslr::NAACproxydata %>% dplyr::filter(Site == "Cedar Island")

## ----loadspt, eval = TRUE-----------------------------------------------------
CedarIslandNC_input <- reslr_load(
  data = CedarIslandNC,
  include_tide_gauge = FALSE,
  include_linear_rate = FALSE,
  TG_minimum_dist_proxy = FALSE,
  list_preferred_TGs = NULL,
  all_TG_1deg = FALSE,
  prediction_grid_res = 50,
  sediment_average_TG = 10
)

## ----dataspt, eval=TRUE-------------------------------------------------------
data <- CedarIslandNC_input$data

## ----datagridspt, eval = TRUE-------------------------------------------------
data_grid <- CedarIslandNC_input$data_grid

## ----printspt, eval=TRUE------------------------------------------------------
print(CedarIslandNC_input)

## ----plotspt,fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(
  x = CedarIslandNC_input,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_proxy_records = TRUE,
  plot_tide_gauges = FALSE
)

## ----runspt,eval = TRUE,message=FALSE,results='hide'--------------------------
res_ni_spline_t <- reslr_mcmc(
  input_data = CedarIslandNC_input,
  model_type = "ni_spline_t",
  CI = 0.95
)

## ----printresspt, eval=TRUE---------------------------------------------------
print(res_ni_spline_t)

## ----summaryspt, eval = TRUE--------------------------------------------------
summary(res_ni_spline_t)

## ---- runsptmore,eval = FALSE-------------------------------------------------
#  res_ni_spline_t <- reslr_mcmc(
#    input_data = CedarIslandNC,
#    model_type = "ni_spline_t",
#    # Update these values
#    n_iterations = 6000, # Number of iterations
#    n_burnin = 1000, # Number of iterations to discard at the beginning
#    n_thin = 4, # Reduces number of output samples to save memory and computation time
#    n_chains = 3 # Number of Markov chains
#  )

## ---- plotsptres,fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(res_ni_spline_t,
  plot_type = "model_fit_plot",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)"
)

## ----plotsptresrate, fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(res_ni_spline_t,
  plot_type = "rate_plot",
  xlab = "Year (CE)",
  y_rate_lab = "Rate of Change (mm per year)"
)

## ----outdfspt, eval = TRUE----------------------------------------------------
output_dataframes <- res_ni_spline_t$output_dataframes
head(output_dataframes)

## ----data2sites, eval = TRUE--------------------------------------------------
# For 2 site
multi_site <- reslr::NAACproxydata %>%
  dplyr::filter(Site %in% c("Cedar Island", "Nassau"))

## ----loadspst, eval = TRUE----------------------------------------------------
multi_site_input <- reslr_load(
  data = multi_site,
  include_tide_gauge = FALSE,
  include_linear_rate = FALSE,
  TG_minimum_dist_proxy = FALSE,
  list_preferred_TGs = NULL,
  all_TG_1deg = FALSE,
  prediction_grid_res = 50,
  sediment_average_TG = 10
)

## ----dataspst,eval=TRUE-------------------------------------------------------
data <- multi_site_input$data
head(data)

## ----datagridspst, eval = TRUE------------------------------------------------
data_grid <- multi_site_input$data_grid
head(data_grid)

## ----printspst, eval=TRUE-----------------------------------------------------
print(multi_site_input)

## ----plotspst,fig.align = 'center',fig.width = 7,fig.height = 10,eval = TRUE----
plot(
  x = multi_site_input,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_proxy_records = TRUE,
  plot_tide_gauges = FALSE
)

## ----runspst,eval = FALSE,results='hide',message=FALSE------------------------
#  res_ni_spline_st <- reslr_mcmc(
#    input_data = multi_site_input,
#    model_type = "ni_spline_st",
#    CI = 0.95
#  )

## ----printspstout,eval=FALSE--------------------------------------------------
#  print(res_ni_spline_st)

## ----summaryspst, eval = FALSE------------------------------------------------
#  summary(res_ni_spline_st)

## ---- runspstmore,eval = FALSE------------------------------------------------
#  res_ni_spline_st <- reslr::reslr_mcmc(
#    input_data = multi_site_input,
#    model_type = "ni_spline_st",
#    # Update these values
#    n_iterations = 6000, # Number of iterations
#    n_burnin = 1000, # Number of iterations to discard at the beginning
#    n_thin = 4, # Reduces number of output samples to save memory and computation time
#    n_chains = 3 # Number of Markov chains
#  )

## ----plotspstres, eval = FALSE,fig.align = 'center',fig.width = 7,fig.height = 5----
#  plot(res_ni_spline_st,
#    plot_type = "model_fit_plot",
#    xlab = "Year (CE)",
#    ylab = "Relative Sea Level (m)"
#  )

## ----plotspstres_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotspstres-1.png"
knitr::include_graphics(url)

## ----plotspstresrate, fig.align = 'center',fig.width = 7,fig.height = 5,eval =FALSE,results='hide',message=FALSE----
#  plot(res_ni_spline_st,
#    plot_type = "rate_plot",
#    xlab = "Year (CE)",
#    y_rate_lab = "Rate of Change (mm per year)"
#  )

## ----plotspstresrate_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotspstresrate-1.png"
knitr::include_graphics(url)

## ----outdfspst, eval = FALSE--------------------------------------------------
#  output_dataframes <- res_ni_spline_st$output_dataframes

## ----data2sitesmore, eval = TRUE----------------------------------------------
# For 9 site
multi_9_sites <- reslr::NAACproxydata %>%
  dplyr::filter(Site %in% c(
    "Cedar Island", "Nassau", "Snipe Key",
    "Placentia", "Cape May Courthouse", "East River Marsh",
    "Fox Hill Marsh", "Swan Key", "Big River Marsh"
  ))

## ----loadnigam, eval = TRUE---------------------------------------------------
multi_9_sites_input <- reslr_load(
  data = multi_9_sites,
  include_tide_gauge = TRUE,
  include_linear_rate = TRUE,
  TG_minimum_dist_proxy = FALSE,
  list_preferred_TGs = NULL,
  all_TG_1deg = TRUE,
  prediction_grid_res = 50,
  sediment_average_TG = 10
)

## ----datanigam,eval=FALSE-----------------------------------------------------
#  data <- multi_9_sites_input$data

## ----datagridnigam, eval = FALSE----------------------------------------------
#  data_grid <- multi_9_sites_input$data_grid

## ----printnigam, eval=TRUE----------------------------------------------------
print(multi_9_sites_input)

## ---- plotnigam,fig.align = 'center',fig.width = 7,fig.height = 10,eval = TRUE----
plot(
  x = multi_9_sites_input,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_proxy_records = TRUE,
  plot_tide_gauges = TRUE
)

## ----runnigam,eval = FALSE,results='hide',message=FALSE-----------------------
#  res_ni_gam_decomp <- reslr_mcmc(
#    input_data = multi_9_sites_input,
#    model_type = "ni_gam_decomp",
#    CI = 0.95
#  )

## ----printnigamout, eval=FALSE------------------------------------------------
#  print(res_ni_gam_decomp)

## ----summarynigam, eval = FALSE-----------------------------------------------
#  summary(res_ni_gam_decomp)

## ---- runnigammore,eval = FALSE-----------------------------------------------
#  res_ni_gam_decomp <- reslr_mcmc(
#    input_data = multi_9_sites_input,
#    model_type = "ni_gam_decomp",
#    # Update these values
#    n_iterations = 6000, # Number of iterations
#    n_burnin = 1000, # Number of iterations to discard at the beginning
#    n_thin = 4, # Reduces number of output samples to save memory and computation time
#    n_chains = 3 # Number of Markov chains
#  )

## ----plotnigamres, eval = FALSE,fig.align = 'center',fig.width = 7,fig.height = 5----
#  plot(res_ni_gam_decomp,
#    plot_type = "model_fit_plot",
#    plot_tide_gauge = FALSE
#  )

## ----plotnigamres_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
knitr::include_graphics("https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotnigamres-1.png")

## ----plotnigamresrate, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_ni_gam_decomp,
#    plot_type = "rate_plot"
#  )

## ----plotnigamresrate_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
knitr::include_graphics("https://raw.githubusercontent.com/maeveupton/reslr/main/docs/reslrvigplots/plotnigamresrate-1.png")

## ----totaldf, eval = FALSE----------------------------------------------------
#  total_model_fit_df <- res_ni_gam_decomp$output_dataframes$total_model_fit_df

## ----plotnigamregres, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_ni_gam_decomp, plot_type = "regional_plot")

## ----plotnigamregres_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
knitr::include_graphics("https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotnigamregres-1.png")

## ----regdf, eval = FALSE------------------------------------------------------
#  regional_component_df <- res_ni_gam_decomp$output_dataframes$regional_component_df

## ----plotnigamregresrate, fig.align = 'center',eval = FALSE-------------------
#  plot(res_ni_gam_decomp, plot_type = "regional_rate_plot")

## ----plotnigamregresrate_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
knitr::include_graphics("https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotnigamregresrate-1.png")

## ----plotnigamlinres, fig.align = 'center',eval = FALSE-----------------------
#  plot(res_ni_gam_decomp, plot_type = "linear_local_plot")

## ----plotnigamlinres_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotnigamlinres-1.png"
knitr::include_graphics(url)

## ----linlocdf, eval = FALSE---------------------------------------------------
#  lin_loc_component_df <- res_ni_gam_decomp$output_dataframes$lin_loc_component_df

## ----linlocrate, eval = FALSE-------------------------------------------------
#  lin_loc_component_rates <- lin_loc_component_df %>%
#    dplyr::group_by(SiteName) %>%
#    dplyr::summarise(
#      linear_rate = unique(linear_rate),
#      linear_rate_err = unique(linear_rate_err)
#    )

## ----plotnigamnonlinres, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_ni_gam_decomp, plot_type = "non_linear_local_plot")

## ----plotnigamnonlinres_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotnigamnonlinres-1.png"
knitr::include_graphics(url)

## ----non_lindf, eval = FALSE--------------------------------------------------
#  non_lin_loc_component_df <- res_ni_gam_decomp$output_dataframes$non_lin_loc_component_df

## ----plotnigamnonlinresrate, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE----
#  plot(res_ni_gam_decomp, plot_type = "non_linear_local_rate_plot")

## ----plotnigamnonlinresrate_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotnigamnonlinresrate-1.png"
knitr::include_graphics(url)

## ----plotnigamall, fig.align = 'center',fig.width = 7,fig.height = 5,eval = FALSE,message=FALSE,results='hide'----
#  plot(res_ni_gam_decomp, plot_type = "nigam_component_plot")

## ----plotnigamall_load, eval = TRUE,echo=FALSE, fig.align = 'center',out.width="100%"----
url <- "https://raw.githubusercontent.com/maeveupton/reslr/main/reslrvigplots/plotnigamall-1.png"
knitr::include_graphics(url)

