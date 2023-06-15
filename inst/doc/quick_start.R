## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"#,
  #fig.path = "qs_vig/"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ---- eval = FALSE, results='hide',message=FALSE------------------------------
#  # Not on CRAN yet
#  #install.packages("reslr")
#  #devtools::install_github("maeveupton/reslr")
#  install_github("maeveupton/reslr")

## ---- start, message=FALSE----------------------------------------------------
library(reslr)

## ---- run_data_qs, eval = TRUE------------------------------------------------
path_to_data <- system.file("extdata", "one_data_site_ex.csv", package = "reslr")
example_one_datasite <- read.csv(path_to_data)

## ---- loadslr_qs, eval = TRUE, message=FALSE, results='hide'------------------
example_one_site_input <- reslr_load(
  data = example_one_datasite)

## ---- plotdata_qs,fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(
  x = example_one_site_input,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_tide_gauges = FALSE,
  plot_caption = TRUE
)

## ---- runslr_qs,eval = TRUE,message=FALSE, results='hide'---------------------
res_one_site_example <- reslr_mcmc(
  input_data = example_one_site_input,
  model_type = "ni_spline_t",
  CI = 0.95
)

## ---- summaryslr_qs, eval = TRUE----------------------------------------------
summary(res_one_site_example)

## ---- plotres_qs, fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(res_one_site_example,
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_type = "model_fit_plot"
)

## ---- plotresrate_qs, fig.align = 'center',fig.width = 7,fig.height = 5,eval = TRUE----
plot(res_one_site_example,
  plot_type = "rate_plot"
)

## ---- dataframeslrres_qs, eval = TRUE-----------------------------------------
output_dataframes <- res_one_site_example$output_dataframes
head(output_dataframes)

