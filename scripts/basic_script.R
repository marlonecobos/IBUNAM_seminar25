################################################################################
# EXAMPLE EXERCISE USING kuenm2 WITH ALGORITHM MAXNET
# Authors: Weverton Trindade, Marlon E. Cobos
# Date updated: 23/02/2025
################################################################################



# Packages ---------------------------------------------------------------------
# Install packages
remotes::install_github("marlonecobos/kuenm2")
install.packages("leaflet")

# Load packages
library(kuenm2)
library(terra)
library(leaflet)
# ------------------------------------------------------------------------------



# Working directory ------------------------------------------------------------
# If needed, set your working directory
#setwd("YOUR/DIRECTORY")
# ------------------------------------------------------------------------------



# Data -------------------------------------------------------------------------
# Dataframe with occurrences (longitude and latitude)
data(occ_data, package = "kuenm2") #Data example
head(occ_data)
str(occ_data)

# Visualize occurrences
pts <- vect(occ_data, geom = c(x = "x", y = "y"), crs = "+init=epsg:4326")
plet(pts, col = "black")

# Raster with environmental variables
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))  # example data

x11()
plot(var)  # notice, these are layers that have been masked to calibration area

# Use only continuous variables
var <- var[[1:4]]
# ------------------------------------------------------------------------------



# Data preparation -------------------------------------------------------------
# Prepare data for maxnet model
sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
                       species = occ_data[1, 1], x = "x", y = "y",
                       raster_variables = var, n_background = 500,
                       features = c("l", "lq", "lqp"),
                       reg_mult = c(0.1, 1))

sp_swd

# Calibration data (Sample with data - SWD)
head(sp_swd$calibration_data)

# Candidate models
head(sp_swd$formula_grid)

# Number of candidate models
nrow(sp_swd$formula_grid)


# Explore the spatial distribution of occurrence and background points
x11()
pbg <- explore_calibration_geo(data = sp_swd, spat_variables = var[[1]],
                               plot = TRUE)

# Explore variable distribution for occurrence and background points
x11()
explore_calibration_hist(data = sp_swd, color_background = "#0000FF80",
                         color_presence = "#FF000080", mfrow = c(2, 2),
                         plot_median = TRUE, breaks = "Scott")
# ------------------------------------------------------------------------------



# Calibrate maxnet models ------------------------------------------------------
m <- calibration(data = sp_swd,  # prepare_data output
                 omission_rate = c(5, 10),  # values (in %) used to calculate the omission rate.
                 omrat_threshold = 10,  # Omission rate used to select the best models
                 parallel = TRUE,  # To run candidate models in parallel
                 ncores = 4,  # Number of cores for parallel processing
                 test_concave = TRUE,  # Whether to test and remove concave curves
                 progress_bar = TRUE)  # Whether to show a progress bar
m

# See summary results for all candidate models
View(m$calibration_results$Summary)

# See results for selected models
View(m$selected_models)

# Select other set of models that perform well among all candidates
# Use another omission rate available in the calibration_results and higher AIC
colnames(m$calibration_results$Summary)

m$omission_rate

m_10 <- sel_best_models(calibration_results = m,  # calibration output
                        algorithm = m$algorithm,
                        omrat_threshold = 5,  # New omission rate
                        delta_aic = 10)  # New Delta AIC

# Compare selected models
View(m_10$selected_models)
# ------------------------------------------------------------------------------



# Model exploration ------------------------------------------------------------
# Fit selected models
fm <- fit_selected(calibration_results = m,  # calibration output
                   n_replicates = 4,  # Number of replicates
                   rep_type = "kfold",  # Replicate type
                   progress_bar = FALSE)
fm


# Variable response curves for best model(s)
## See variables available to explore response curves
names(fm$Models$Model_35$Full_model$betas)

## Plots
x11()
par(mfrow = c(1, 3), cex = 0.8)
response_curve(models = fm, variable = "bio_1")
response_curve(models = fm, variable = "bio_7")
response_curve(models = fm, variable = "bio_12")


# Two-Way interaction response curves
## Plot with our example
x11()
resp2var(models = fm, modelID = "Model_35",
         variable1 = "bio_1", variable2 = "bio_12")

## An example included in the package of a model with product responses
data("fitted_model_glmnet", package = "kuenm2")

## Response curves (same variables but model includes product)
x11()
resp2var(models = fitted_model_glmnet, modelID = "Model_13",
         variable1 = "bio_1", variable2 = "bio_12")


# Variable importance (variable contribution to model)
## Asses importance
imp_maxnet <- var_importance(models = fm)

## Plot importance (derived directly from package enmpa)
x11()
plot_importance(imp_maxnet)
# ------------------------------------------------------------------------------



# Model transfers --------------------------------------------------------------
# Predict selected models for a single scenario
## Prediction in M
p <- predict_selected(models = fm,  # fit_selected output
                      raster_variables = var,  # SpatRaster of variables for predicting
                      consensus_per_model = TRUE,  # compute consensus for each model across its replicates
                      consensus_general = TRUE,  # compute a general consensus across all models
                      consensus = c("median", "range", "mean", "stdev"),
                      clamping = FALSE, type = "cloglog", progress_bar = FALSE)

# Plot replicates
x11()
plot(p$Model_35$Replicates)

# Plot consensus between replicates
x11()
plot(p$Model_35$Model_consensus)

# Plot consensus between models
x11()
plot(p$General_consensus)


# Predictions to multiple scenarios
# Organize and structure future climate variables from WorldClim ####
# Set the input directory containing the raw future climate variables.
# For this example, the data provided with the package will be used.
in_dir <- system.file("extdata", package = "kuenm2")

list.files(in_dir)  # Raw variables downloaded from WorldClim

# Create a "_variables" folder in your directory to save results
out_dir_future <- file.path("Future_variables")
dir.create(out_dir_future)

# Use this function to Organize and rename the future climate data, structuring it by year and GCM
organize_future_worldclim(input_dir = in_dir, output_dir = out_dir_future,
                          name_format = "bio_", variables = NULL,
                          fixed_variables = NULL, mask = NULL, overwrite = TRUE)
# Check the folder


# Preparation of data for model projections
# We already have the variables of the future organized with organize_future_worldclim
# We also need the raw variables of the present in a folder
# Let's save the present variables (var) in a folder
out_dir_current <- "Current_raw"
dir.create(out_dir_current)

# Save current variables in this directory
writeRaster(var, file.path(out_dir_current, "Variables.tif"))


# Now, we can prepare projections data using fitted models to check variables
pr <- prepare_proj(models = fm, present_dir = out_dir_current, past_dir = NULL,
                   past_period = NULL, past_gcm = NULL,
                   future_dir = out_dir_future,
                   future_period = c("2041-2060", "2081-2100"),
                   future_pscen = c("ssp126", "ssp585"),
                   future_gcm = c("ACCESS-CM2", "MIROC6"),
                   write_file = FALSE, filename = NULL,
                   raster_pattern = ".tif*")
pr
str(pr)

# Create folder to save projection results
out_dir <- file.path("Projection_results/maxnet")
dir.create(out_dir, recursive = TRUE)


# Project selected models for multiple scenarios
p <- project_selected(models = fm, projection_data = pr, out_dir = out_dir,
                      consensus_per_model = TRUE, consensus_general = TRUE,
                      consensus = c("median", "range", "mean", "stdev"),
                      write_replicates = TRUE, clamping = FALSE,
                      var_to_clamp = NULL, type = "cloglog", overwrite = TRUE,
                      parallel = FALSE, ncores = 1, parallelType = "doSNOW",
                      progress_bar = TRUE, verbose = TRUE)

View(p$paths)

## Files created
list.files(out_dir, recursive = TRUE)


## Import some projection
res_present <- rast("Projection_results/maxnet/Present/Present/General_consensus.tiff")
res_future <- rast("Projection_results/maxnet/Future/2081-2100/ssp585/MIROC6/General_consensus.tiff")

## Plot projections
x11()
par(mfrow = c(1, 2))
plot(res_present$mean, main = "Present")
plot(res_future$mean, main = "2081-2100 - SSP85 - MIROC 6")
# ------------------------------------------------------------------------------



# Post-modeling analyses -------------------------------------------------------
# Compute areas of contraction, expansion and stability in the future
changes <- proj_changes(model_projections = p, reference_id = 1,
                        consensus = "mean", include_id = NULL, user_thr = NULL,
                        by_gcm = TRUE, by_change = TRUE, general_summary = TRUE,
                        force_resample = TRUE, write_results = FALSE,
                        output_dir = NULL, overwrite = FALSE,
                        write_bin_models = FALSE, return_rasters = TRUE)

## Plot results
x11()
plot(changes$Binarized)  # SpatRaster with the binarized models for each GCM
plot(changes$Results_by_gcm)  # SpatRaster with the computed changes for each GCM
plot(changes$Results_by_change$`Future_2041-2060_ssp126`)  # List containing the SpatRaster with each computed change for each GCM
plot(changes$Summary_changes)  # SpatRaster with the general summary


# Variance in predictions coming from distinct sources in ENMs
v <- modvar(model_projections = p, by_replicate = T, by_gcm = TRUE,
            by_model = TRUE, consensus = "median", write_files = FALSE,
            output_dir = NULL, return_rasters = TRUE, progress_bar = FALSE,
            verbose = TRUE, overwrite = FALSE)

x11()
plot(v$Present$by_rep)  # Variance coming from replicates in Present projection
plot(v$Present$by_model)  # Variance coming from models in Present projection
plot(v$`Future_2041-2060_ssp126`$by_rep)  # Variance coming from replicates in one of the future projections
plot(v$`Future_2041-2060_ssp126`$by_model)  # Variance coming from models in one of the future projections

plot(v$`Future_2041-2060_ssp126`$by_gcm)  # Variance coming from GCMs in one of the future projections
plot(v$`Future_2041-2060_ssp585`$by_gcm)  # Variance coming from GCMs in one of the future projections
plot(v$`Future_2081-2100_ssp126`$by_gcm)  # Variance coming from GCMs in one of the future projections
plot(v$`Future_2081-2100_ssp585`$by_gcm)  # Variance coming from GCMs in one of the future projections


# Analysis of extrapolation risks using the MOP metric
## Create folder to save MOP results
out_dir <- file.path("MOP_results")
dir.create(out_dir)

## Run MOP
kmop <- kuenm_mop(data = sp_swd, subset_variables = TRUE, mask = NULL,
                  fitted_models = fm, projection_data = pr, out_dir = out_dir,
                  type = "detailed", calculate_distance = FALSE,
                  where_distance = "in_range", distance = "euclidean",
                  scale = FALSE, center = FALSE, fix_NA = TRUE, percentage = 1,
                  comp_each = 2000, tol = NULL, rescale_distance = FALSE,
                  parallel = FALSE, n_cores = 1, progress_bar = FALSE,
                  overwrite = TRUE)
kmop

## Import simple mop
simple_mop_files <- list.files("MOP_results/Future/",
                               recursive = TRUE, full.names = TRUE,
                               pattern = "simple")

## Remove aux files
simple_mop_files <- simple_mop_files[!grepl("aux", simple_mop_files)]

## Read rasters
simple_mop <- rast(simple_mop_files)

## Rename rasters
names(simple_mop) <- sub("\\_mopsimple.tif$", "", gsub("/", "-",
                                              sub(".*MOP_results/", "",
                                                  simple_mop_files)))

## Plot MOP results
x11()
plot(simple_mop)


## Import detailed mop
detailed_mop_files <- list.files("MOP_results/Future/",
                                 recursive = TRUE, full.names = TRUE,
                                 pattern = "combined")

## Remove aux and csv files
detailed_mop_files <- detailed_mop_files[!grepl("aux|csv", detailed_mop_files)]

## Read rasters
detailed_mop <- rast(detailed_mop_files)

## Rename rasters
names(detailed_mop) <- sub("\\.tif$", "", gsub("/", "-",
                                               sub(".*MOP_results/", "",
                                                   detailed_mop_files)))
## Remove results with NA
na_results <- sapply(detailed_mop, function(x) !is.na(minmax(x)[1]))
detailed_mop <- detailed_mop[[na_results]]


## Plot first 8 mop results
x11()
plot(detailed_mop, cex.main= 0.5)
# ------------------------------------------------------------------------------
