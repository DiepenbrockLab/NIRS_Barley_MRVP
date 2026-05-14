#Title: Custom NIRS calibration for Barley F2 populations

#Clear the environment and load packages

#{r call libraries}
rm(list=ls())
library(prospectr)
library(pls)
library(ggplot2)
library(ggfortify)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggrepel)
library(magrittr)
library(robustbase)
library(glue)
library(R.utils)
set.seed(12345)

#{r exclude outliers?}
# Exclude spectral and wet chem outliers?
Spec_exclude_outliers = TRUE
WetChem_exclude_outliers = TRUE

#{r read macronutrients}

# Read FOSS Cereals/Grains predicted macronutrients 

setwd("C:/Users/chdiep/Documents/ResearchProjects/BARD/AbelinaCode_20260106")
BarleyF2.whole.macs = read.csv("F2_whole_macros_cleaned_updated.csv",head=T)

df <- BarleyF2.whole.macs %>%
  mutate(location = case_when(
    startsWith(Sample.Number, "Y_") ~ "Yotvata",
    startsWith(Sample.Number, "M_") ~ "Mibhor",
    TRUE ~ "Unknown"  # default value for unexpected prefixes
  ))
BarleyF2.whole.macs$Sample.Type = "Whole"

#{r read and format spectra}

########
# Read spectra (Spectragryph)
BarleyF2.whole.spectra = read.csv("F2_whole_scanned_spectra_cleaned.csv")

tail(colnames(BarleyF2.whole.spectra))
if(colnames(BarleyF2.whole.spectra)[ncol(BarleyF2.whole.spectra)] == "X.1"){
  BarleyF2.whole.spectra = BarleyF2.whole.spectra[,-ncol(BarleyF2.whole.spectra)]
}

# Remove extra row from spectra, transpose
F2.spectra = t(BarleyF2.whole.spectra[-1,])

# Improve column and row names and drop first column
No.decimal = gsub("\\.5.*", ".5", F2.spectra[1,])
No.decimal = gsub("\\.0.*", ".0", No.decimal)
colnames(F2.spectra) = paste0("Band_", No.decimal)
F2.spectra.cln = F2.spectra[-1,]

# Check sample names to make sure they are correctly aligned
rownames(F2.spectra.cln)[1:100] # does this one match with
BarleyF2.whole.macs$Sample.Number[1:100] # this one?

#{r extract plot numbers from sample names}

# **Check to make sure  sample names are correctly aligned in prev. chunk**
# Repair truncated names (this steps will vary based on sample set)
rownames(F2.spectra.cln) <- sub(".*(M_|Y_)", "\\1", rownames(F2.spectra.cln))

#Setting rownames to match sample.number
rownames(BarleyF2.whole.macs) <- BarleyF2.whole.macs$Sample.Number

#Check structure
str(F2.spectra.cln)

#Convert to numerical matrix without altering dimensions
storage.mode(F2.spectra.cln) <- "numeric"
str(F2.spectra.cln)

dat_spectral = F2.spectra.cln

#{r read and format wetchem}
#Read in wetchem data and format dataframe
list.files()
dat_wetchem = read.csv("F2_anlab_results_without_dup.csv")
str(dat_wetchem)

#Remove extra columns
dat_wetchem <- dat_wetchem[, !names(dat_wetchem) %in% "sample"]
dat_wetchem$Moisture = 100 - dat_wetchem$DM
colnames(dat_wetchem)

names(dat_wetchem)[3:7] = capitalize(names(dat_wetchem)[3:7])
colnames(dat_wetchem)

#Remove DM column
dat_wetchem = dat_wetchem[,c(1,3:7)]

#{r calculate spectral pretreatments}

# Pull out wavelengths / lambdas
lambdas <- as.numeric(gsub("Band_", "", colnames(dat_spectral)))

# Set window size for SG pretreatments and create trimmed lambdas
w <- 11
#trim <- (w - 1) / 2
#lam_trim <- lambdas[(trim + 1):(length(lambdas) - trim)] #was not being used for anything

# Calculate spectral Pretreatments
# Standard Normal Variate (SNV)
snv <- standardNormalVariate(dat_spectral)

# Multiplicative Scatter Correction (MSC)
msc <- msc(dat_spectral)

# Detrending (DT)
dt <- detrend(dat_spectral, wav = lambdas)

# Baseline Correction (BL)
bl <- baseline(dat_spectral, wav = lambdas)
colnames(bl) = paste0("Band_",format(as.numeric(colnames(bl)), nsmall = 1))

# Apply the Savitzky-Golay filter for smoothing (0th derivative)
sg <- savitzkyGolay(dat_spectral, p = 3, w = w, m = 0)#  polynomial degrees (p, 2 and 3 seems to be most common), windowsize (w), order of differentiation (m)

# Compute the first derivative using the Savitzky-Golay filter
sg1 <- savitzkyGolay(dat_spectral, p = 3, w = w, m = 1)

# Compute the second derivative using the Savitzky-Golay filter
sg2 <- savitzkyGolay(dat_spectral, p = 3, w = w, m = 2)

# Compute SNV of detrending filter
snv.dt <- detrend(snv, wav = lambdas)

# Compute the first derivative using the Savitzky-Golay filter
snv.sg <- savitzkyGolay(snv, p = 3, w = w, m = 0)

# Compute the first derivative using the Savitzky-Golay filter
snv.sg1 <- savitzkyGolay(snv, p = 3, w = w, m = 1)

# Compute the second derivative using the Savitzky-Golay filter
snv.sg2 <- savitzkyGolay(snv, p = 3, w = w, m = 2)

# p is the polynomial order, w is the window size (must be odd), and m is the derivative order

# # Inspect the results of pretreatments
matplot(lambdas, t(dat_spectral), type = "l", lty = 1, col = 1:nrow(dat_spectral),xlab = "Wavelength", ylab = "Absorbance; log(1/R)", main = "No Pretreatment")
matplot(lambdas, t(snv), type = "l", lty = 1, col = 1:nrow(snv),xlab = "Wavelength", ylab = "Transformed Absorbance", main = "Standard Normal Variate (SNV)")
matplot(lambdas, t(msc), type = "l", lty = 1, col = 1:nrow(msc),xlab = "Wavelength", ylab = "Transformed Absorbance", main = "Multiplicative Scatter Correction (MSC)")
matplot(lambdas, t(dt), type = "l", lty = 1, col = 1:nrow(dt),xlab = "Wavelength", ylab = "Transformed Absorbance", main = "Detrending (DT)")
matplot(lambdas, t(bl), type = "l", lty = 1, col = 1:nrow(bl),xlab = "Wavelength", ylab = "Transformed Absorbance", main = "Baseline (BL)")
matplot(lambdas, t(snv.dt), type = "l", lty = 1, col = 1:nrow(snv.dt),xlab = "Wavelength", ylab = "Transformed Absorbance", main = "Detrended Standard Normal Variate (SNV.DT)")

#Make list of spectral pretreatments and untreated spectra
pretreats <- list(spec = dat_spectral, snv = snv, msc = msc,  bl = bl, dt= dt, snv.dt = snv.dt, sg = sg, sg1 = sg1, sg2 = sg2,  snv.sg = snv.sg, snv.sg1 = snv.sg1, snv.sg2 = snv.sg2) # Note: SG and SNV-SG variants have 10 fewer wavelengths due to windowing (w = 11)

#{r compare number of columns across pretreatments}
#Number of columns differ across pretreatments, this chunk checks how many and which ones
lapply(pretreats, ncol)

# Function to compare column names with the first data frame
compare_colnames <- function(df, reference_df) {
  list(
    only_in_reference = setdiff(colnames(reference_df), colnames(df)),
    only_in_current = setdiff(colnames(df), colnames(reference_df))
  )
}

# Apply the function to each data frame in the list
# Compare each data frame to the first one in the list, which is the raw spectra
results <- lapply(pretreats[-1], compare_colnames, reference_df = pretreats[[1]])

#{r PCA on pretreatments}

#Define PCA function
PCA <- function(X_list, name) {
  cat("\nProcessing pretreatment:", name, "\n")
  
  # Remove constant columns
  constant_cols <- which(apply(X_list, 2, var) < 1e-10)
  if (length(constant_cols) > 0) {
    warning("Percent of columns removed for low variance: ", 
            round(length(constant_cols) / ncol(X_list) * 100, 2), "%")
    X_list <- X_list[, -constant_cols, drop = FALSE]
  }
  
  # Check if the dataset is now empty
  if (ncol(X_list) == 0) {
    warning("All columns were constant in ", name, " - Skipping PCA.")
    return(NULL)
  }
  
  # Perform PCA
  pc_res <- tryCatch(
    prcomp(X_list, scale. = TRUE),
    error = function(e) {
      warning("PCA failed for ", name, " with error: ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(pc_res)) {
    # Create environment labels
    env_labels <- substr(rownames(X_list), 1, 5)
    
    # PCA scatter plot
    print(
      autoplot(pc_res, label = TRUE, label.size = 2, main = name, geom = "text",
               data = data.frame(env = env_labels), colour = 'env')
    )
    
    # Scree plot
    var_explained <- pc_res$sdev^2 / sum(pc_res$sdev^2)
    scree_df <- data.frame(PC = seq_along(var_explained), Variance = var_explained)
    
    print(
      ggplot(scree_df, aes(x = PC, y = Variance)) +
        geom_line() + geom_point() +
        scale_y_continuous(labels = scales::percent) +
        labs(title = paste("Scree Plot -", name), x = "Principal Component", y = "Variance Explained") +
        theme_minimal()
    )
  }
  return(pc_res)
}

# Apply the PCA function and store results in a list
pca_results <- lapply(seq_along(pretreats), function(i) 
  PCA(pretreats[[i]], names(pretreats)[i])
)

#{r Outlier detection spectral pretreatments}

#Define function to identify and remove outliers in spectral pretreatments
outlier_spectral <- function(pretreatments, pca_results, removal_threshold = round(length(pretreats)*0.75)) { 
  
  outlier_results <- list()  # To store outlier info for each pretreatment
  
  for (i in seq_along(pretreatments)) {
    name <- names(pretreatments)[i]
    
    # Get PCA scores for the current pretreatment (using first 5 PCs)
    pc_scores <- pca_results[[i]]$x[, 1:5]
    
    # Calculate robust Mahalanobis distances
    mcd <- covMcd(pc_scores)
    mahal_distances <- mahalanobis(pc_scores, mcd$center, mcd$cov)
    
    # Use a fixed threshold at the 99.5th percentile
    perc_thresh <- quantile(mahal_distances, 0.995)
    
    # Identify outliers (using rownames as sample IDs)
    outliers <- rownames(pc_scores)[mahal_distances > perc_thresh]
    n <- nrow(pc_scores)
    
    # Store outlier info for this pretreatment
    outlier_results[[name]] <- list(
      outliers = outliers,
      mahalanobis = mahal_distances,
      threshold = perc_thresh,
      variance_explained = summary(pca_results[[i]])$importance[2, 1:5]
    )
  }
  
  # Combine outlier results across pretreatments
  all_outliers <- unlist(lapply(outlier_results, function(x) x$outliers))
  outlier_counts <- table(all_outliers)
  cat("\nOutlier frequency across pretreatments:\n")
  print(outlier_counts)
  
  # Determine which samples are flagged in removal_threshold or more pretreatments
  samples_to_remove <- names(outlier_counts)[outlier_counts >= removal_threshold]
  cat("\nSamples to remove (flagged in", removal_threshold, "or more pretreatments):\n")
  print(samples_to_remove)
  
  # If removal is enabled, remove these samples from every pretreatment dataset
  if(Spec_exclude_outliers){
    cleaned_data <- lapply(pretreatments, function(dt) {
      dt[!rownames(dt) %in% samples_to_remove, ]
    })
  } else {
    cleaned_data <- pretreatments
    cat("Outliers retained for all pretreatments.\n")
  }
  
  return(list(
    cleaned_data = cleaned_data,  # Pretreatment data with joint outlier removal
    outlier_info = outlier_results,  # Outlier info per pretreatment
    removal_list = samples_to_remove # List of samples removed
  ))
}

# Run the function using both pretreatments and PCA results
results <- outlier_spectral(pretreats, pca_results)

# Access cleaned data
cleaned_data <- results$cleaned_data

# View outliers for a specific pretreatment (e.g., "SG.2")
all_pretre_outs <- unlist(lapply(results$outlier_info, function(x) x$outliers))

#{r view spectral outliers in biplot}
# List of samples flagged for removal across pretreatments
samples_to_remove <- results$removal_list

# Loop through each pretreatment
for (i in seq_along(results$outlier_info)) {
  pretreatment_name <- names(results$outlier_info)[i]
  print(pretreatment_name)
  
  # Get PCA results for the current pretreatment
  pca_result <- pca_results[[i]]
  pca_scores <- pca_result$x[, 1:2]
  
  # Compute variance explained
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  pc1_var <- round(variance_explained[1], 2)
  pc2_var <- round(variance_explained[2], 2)
  
  # Get outliers for the current pretreatment
  outliers_this_pretreat <- results$outlier_info[[pretreatment_name]]$outliers
  
  # Create data frame for plotting
  pca_df <- data.frame(
    PC1 = pca_scores[, 1],
    PC2 = pca_scores[, 2],
    Sample = rownames(pca_scores),
    OutlierPretreat = rownames(pca_scores) %in% outliers_this_pretreat,
    OutlierGlobal = rownames(pca_scores) %in% samples_to_remove
  )
  
  # Create a combined Outlier Status
  pca_df <- pca_df %>%
    mutate(Status = case_when(
      OutlierGlobal ~ "Global Outlier",
      OutlierPretreat ~ "Pretreatment Outlier",
      TRUE ~ "Normal"
    ))
  
  # Plot
  plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Status)) +
    geom_point(size = 3, alpha = .7) +
    scale_color_manual(
      values = c(
        "Pretreatment Outlier" = "orange",
        "Global Outlier" = "red"
      )
    ) +
    labs(
      title = paste("PCA Biplot for", pretreatment_name),
      x = paste("Principal Component 1 (", pc1_var, "%)", sep = ""),
      y = paste("Principal Component 2 (", pc2_var, "%)", sep = ""),
      color = "Outlier Type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(plot)
}

#{r PCA on wet chem data}

dat_wetchem_noMissing = dat_wetchem

pc_res_wetchem = prcomp(dat_wetchem_noMissing[,-1], center=T, scale.=T)
plot(pc_res_wetchem)
autoplot(pc_res_wetchem,label=T,label.size=3.00,loadings=T,loadings.label=T,loadings.label.size=3.75,scale=T)

#{r outlier detection for wetchem}
# Outlier detection for wetchem:
# Function to detect outliers using IQR
detect_outliers <- function(x) {
  outliers <- rep(FALSE, length(x)) # Initialize a vector of FALSEs
  x_non_na <- na.omit(x) # Remove missing values
  Q1 <- quantile(x_non_na, 0.25)
  Q3 <- quantile(x_non_na, 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  outliers[!is.na(x) & (x < lower_bound | x > upper_bound)] <- TRUE # Mark outliers in the non-NA positions
  return(outliers)
}

# Apply the function to each numeric column and add outlier column
WC_outliers <- dat_wetchem %>%
  select(-plot) %>%
  mutate(across(everything(), detect_outliers, .names = "outlier_{.col}"))

# Combine the original data frame with the outlier information
dat_wetchem_outs <- bind_cols(SampleNumber = dat_wetchem$plot, WC_outliers)

# Gather data for visualization
df_long_WC <- dat_wetchem_outs %>%
  pivot_longer(cols = -c(SampleNumber, starts_with("outlier")), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = gsub("_percent", "", Metric)) %>%
  select(-starts_with("outlier"))

# Pivot dat_wetchem_outs for outliers
df_long_outs <- dat_wetchem_outs %>%
  pivot_longer(cols = starts_with("outlier"), names_to = "Metric", values_to = "Outlier") %>%
  mutate(Metric = gsub("^outlier_|_percent", "", Metric)) %>%
  select(SampleNumber, Metric, Outlier)

df_long <- merge(df_long_WC, df_long_outs, by=c("SampleNumber", "Metric"))

# Filter outliers
outliers <- df_long[df_long$Outlier == TRUE, ]

# Plotting
ggplot(na.omit(df_long), aes(y = Value)) +
  geom_boxplot(aes(x = ""), outlier.shape = NA) +  # x = "" ensures one box per facet
  geom_jitter(aes(x = ""), alpha = 0.3) +
  geom_point(data = outliers, aes(x = "", color = "red"), size = 3) +  # Set x = "" for outliers
  geom_text_repel(data = outliers, aes(x = "", label = SampleNumber),
                  box.padding = .5, direction = "y", segment.color = NA,
                  force = 0.5, size = 2.5, color = "red") +
  scale_color_identity() +
  theme_minimal() +
  labs(title = "Wet Chem Outlier Detection",
       y = "Percent Composition (%)",
       color = "Outlier") +
  facet_wrap(~Metric, scales = "free")

#outliers
out.pct.ash = sum(dat_wetchem_outs$outlier_Ash)/nrow(dat_wetchem_outs)*100
out.pct.prot = sum(dat_wetchem_outs$outlier_Protein)/nrow(dat_wetchem_outs)*100
out.pct.moi = sum(dat_wetchem_outs$outlier_Moisture)/nrow(dat_wetchem_outs)*100
out.pct.fat = sum(dat_wetchem_outs$outlier_Fat)/nrow(dat_wetchem_outs)*100

cat("Outliers for Ash are", round(out.pct.ash,1), "% of data.\n",
    "Outliers for Protein are", round(out.pct.prot,1), "% of data.\n",
    "Outliers for Moisture are", round(out.pct.moi,1), "% of data.\n",
    "Outliers for Fat are", round(out.pct.fat,1), "% of data.\n")

#{r outlier removal wetchem data}
#Create new object, copy of wetchem data
dat_wetchem_NoOuts = dat_wetchem
rownames(dat_wetchem_NoOuts) <- dat_wetchem_NoOuts$plot

 #Loop through outliers and replace specified cells with NA
for(i in 1: nrow(outliers)) {
 Trait_i = outliers$Metric[i]
  Sample_i = as.character(outliers$SampleNumber[i])
  
  # Ensure that Sample_i and Trait_i exist in the dataframe
  if (Sample_i %in% rownames(dat_wetchem_NoOuts) && Trait_i %in% colnames(dat_wetchem_NoOuts)) {
    dat_wetchem_NoOuts[Sample_i, Trait_i] <- NA
  }else{
    print("Missing outlier in wetchem.")
  }
}

dat_wetchem_NoOuts
dim(dat_wetchem_NoOuts)
dim(dat_wetchem_noMissing)

#{r merge wetchem and spectral data}
 wetchem_data = dat_wetchem_NoOuts

# Define a function to merge wet-chemistry data with spectral data for a specific pretreatment

merge_wetchem_with_spectra <- function(wetchem_data, spectra_data, pretreatment_name) {
  
 if (WetChem_exclude_outliers == TRUE & Spec_exclude_outliers == TRUE) {
    merged_data <- merge(wetchem_data, spectra_data, by = 'row.names')
  } else if (WetChem_exclude_outliers == TRUE & Spec_exclude_outliers == FALSE) {
    merged_data <- merge(wetchem_data, pretreats[[pretreatment_name]], by = 'row.names')
  } else if (WetChem_exclude_outliers == FALSE & Spec_exclude_outliers == TRUE) {
    merged_data <- merge(wetchem_data, spectra_data, by = 'row.names')
  } else if (WetChem_exclude_outliers == FALSE & Spec_exclude_outliers == FALSE){
    merged_data <- merge(wetchem_data, pretreats[[pretreatment_name]], by = 'row.names')
  } else{
    print(paste0("Error: WetChem_exclude_outliers and Spec_exclude_outliers should
               be logical (TRUE or FALSE) but are:"),
         Spec_exclude_outliers, WetChem_exclude_outliers)
    break
  }
  
   #Set row names and clean up the merged dataset
  rownames(merged_data) <- merged_data$Row.names
  merged_data <- merged_data[, -c(1,2)]  # Adjust this as needed for your final dataset
  
 return(merged_data)
}

# Iterate over all pretreatments in cleaned_data
pretreatment_results <- list()

for (pretreatment_name in names(cleaned_data)) {
  spectra_data <- as.data.frame(cleaned_data[[pretreatment_name]])
  
  merged_result <- merge_wetchem_with_spectra(wetchem_data, spectra_data, pretreatment_name)
  
  pretreatment_results[[pretreatment_name]] <- merged_result
  
  # Display dimensions of merged dataset for sanity check
  cat("Pretreatment", pretreatment_name, "\n Rows:", dim(merged_result)[1], "\n Cols:", dim(merged_result)[2], "\n")
}

# Loop through each merged dataset in pretreatment_results
for (name in names(pretreatment_results)) {
  
  # Extract the merged dataset
  merged_data <- pretreatment_results[[name]]
  
  # Identify samples in wetchem that are not in the spectral dataset
  missing_samples <- setdiff(rownames(wetchem_data), rownames(merged_data))
  
  # Print results for the current pretreatment
  cat("\n--- Pretreatment:", name, "---\n")
  cat("Samples (with wetchem) removed from analysis:\n")
  print(missing_samples)
}

# kenStone produces a list with row index of the points selected for calibration
# Note: the line below is currently using PCs on the spectra rather than the spectra themselves (as the script otherwise wanted # indiv to be > # wavelengths)
# Using 2 PCs as % var explained tapers off quickly (see results of line 23)
# Currently selecting k = 80% of samples to be in the calibration set, leaving 20% in the test set

# Initialize an empty list to store kenStone results for each pretreatment
ken_mahal_results <- list()

for (pretreatment_name in names(pretreatment_results)) {
  
  merged_result <- pretreatment_results[[pretreatment_name]]
  
  # Select numeric columns
  numeric_data <- merged_result[, sapply(merged_result, is.numeric)]
  numeric_data <- as.data.frame(lapply(numeric_data, as.numeric))
  
  # Remove rows with NA, NaN, Inf
  numeric_data <- numeric_data[apply(numeric_data, 1, function(x) all(is.finite(x))), ]
  
  n_samples <- nrow(numeric_data)
  k <- round(0.80 * n_samples)
  
  ken_mahal <- kenStone(numeric_data,
                        k = k,
                        metric = "mahal",
                        pc = 2)
  
  ken_mahal_results[[pretreatment_name]] <- ken_mahal
  
  # Plot
  plot(ken_mahal$pc[, 1],
       ken_mahal$pc[, 2],
       col = rgb(0, 0, 0, 0.3),
       pch = 19,
       xlab = "PC1",
       ylab = "PC2",
       main = paste("Kennard-Stone -", pretreatment_name))
  grid()
  
  points(ken_mahal$pc[ken_mahal$model, 1],
         ken_mahal$pc[ken_mahal$model, 2],
         pch = 19, col = "red")
}

 # If file is already made, run commented lines below. If not, run the function at line 569.
 # ncomp_overrides_df <- read.csv("Ncomp_overrides_2026_01_05_updated.csv", row.names=1) 
 # ncomp_overrides <- lapply(seq_len(ncol(ncomp_overrides_df)), function(j) {
 #   vals <- as.numeric(ncomp_overrides_df[[j]])
 #   names(vals) <- rownames(ncomp_overrides_df)
 #   vals
 # })
 # names(ncomp_overrides) <- colnames(ncomp_overrides_df)

 #Function to select number of components interactively
 selectNcompInteractive <- function(pretreatment_results, traits = 1:5) {
   ncomp_selection <- list()  # initialize the list
   
   # Open external window for plotting (platform-specific)
   if (Sys.info()["sysname"] == "Windows") {
     x11()
   } else if (Sys.info()["sysname"] == "Darwin") {
     quartz()
   } else {
     x11()
   }
   
   for (pt in names(pretreatment_results)) {
   cat("Processing pretreatment:", pt, "\n")
   dt <- as.data.frame(pretreatment_results[[pt]])  # ensure it's a data.frame
   
   # Create a vector to store final number of components for each trait
   ncomp_vector <- numeric(length(traits))
   names(ncomp_vector) <- colnames(dt)[traits]
   
   for (i in traits) {
   trait_name <- colnames(dt)[i]
   cat("  Processing trait:", trait_name, "\n")
   
   # Fit baseline PLSR model using 5-fold CV
   mod <- plsr(as.formula(paste(trait_name, "~ .")),
               data = dt,
               validation = "CV",
               segments = 5)
   
   # Extract RMSEP values
   rmsep_vals <- RMSEP(mod)$val[1, 1, ]
   comps <- 1:length(rmsep_vals)
   
   # Automatic selection
   default_ncomp <- as.numeric(selectNcomp(mod, method = "onesigma"))
   if (default_ncomp == 0) default_ncomp <- 1
   
   # Plot RMSEP vs components
   plot(comps, rmsep_vals, type = "b", pch = 16, col = "blue",
        xlab = "Number of Components", ylab = "RMSEP",
        main = paste("Screeplot for", trait_name, ":", pt))
   text(comps, rmsep_vals, labels = comps, pos = 3, cex = 0.8, col = "gray20")
   points(default_ncomp, rmsep_vals[default_ncomp], col = "red", pch = 19, cex = 1.5)
   abline(v = default_ncomp, col = "red", lty = 2)
   grid()
   flush.console()
   
   cat("Automatic selection using onesigma method...\n")
   cat("Automatically selected number of components for", trait_name, "=", default_ncomp, "\n")
   
   # Prompt override
   user_input <- readline(prompt = "    Enter number of components to override (or press Enter to use default): ")
   final_ncomp <- ifelse(nchar(trimws(user_input)) > 0, as.numeric(user_input), default_ncomp)
   cat("    Final number of components for", trait_name, "=", final_ncomp, "\n\n")
   
   ncomp_vector[i] <- final_ncomp
 }
 
 # Save vector for this pretreatment
 ncomp_selection[[pt]] <- ncomp_vector
 }

return(ncomp_selection)
 }
 
 
 #{r Main Execution Block}
 ext_HD <- "Temp_PLSR_customCal_barley_Model_outputs" 
 setwd(ext_HD)
 
 # Create output directory for today's modeling runs
 dir <- file.path("saved_models", Sys.Date()) 
 dir.create(dir, showWarnings = FALSE, recursive = TRUE)
 
 niter = 5

 for(q in 4:niter){
   
# Example usage (assuming 'pretreats' is your list of pretreatment datasets):
ncomp_overrides <- selectNcompInteractive(pretreatment_results, traits = 1:5)# change to 1:8 once adding starch
print(ncomp_overrides)

# If you made a mistake you can correct like so:
#ncomp_overrides$sg2["Protein"] <- 6 #protein
#ncomp_overrides$snv.sg1[["Moisture"]] <- 19 #moisture
#ncomp_overrides$snv.sg1["Protein"] <- 9 #protein
#ncomp_overrides$sg2 <- ncomp_overrides$sg2[1:5]
# ncomp_overrides$snv.dt[[10]] <- 12 #fat
#print(ncomp_overrides) #check corrections worked

# Write to file for option of skipping this step
write.csv(ncomp_overrides, paste0("Ncomp_overrides_",Sys.Date(),"_updated_iter",q,".csv"))

# PLSR function to apply across pretreatments (no LOO CV)
FitAndSaveModels <- function(pretreatment_results, ken_mahal_results, ncomp_overrides, traits = 1:5, output_dir = dir) {#change to 1:8 when adding starch
  
  # Identify models already saved (CV and KS both exist)
  cv_done <- list.files(output_dir, pattern = "^mod_cv_.*\\.rds$") |> 
    sub("^mod_cv_", "", x = _) |> sub("\\.rds$", "", x = _)
  
  ks_done <- list.files(output_dir, pattern = "^mod_ks_.*\\.rds$") |> 
    sub("^mod_ks_", "", x = _) |> sub("\\.rds$", "", x = _)
  
  completed <- intersect(cv_done, ks_done)
  
  # Loop through each pretreatment
  for (pt in names(pretreatment_results)) {
    cat("Processing Pretreatment:", pt, "\n")
    dt <- as.data.frame(pretreatment_results[[pt]])
    ken_mahal <- ken_mahal_results[[pt]]
    
    # Loop through each trait (within pretreatment)
    for (t in traits) {
      
      trait_name <- colnames(dt)[t]
      task_id <- paste0(pt, "_", gsub("\\s+", "_", trait_name,"_iter",q))
      if (task_id %in% completed) {
        cat(sprintf("  Skipping %s: already completed\n", task_id))
        next
      }
      
      cat("  Processing Trait:", trait_name, "\n")
      
      # Subset and filter rows with complete data for this trait
      # Get the trait name and all predictor columns that start with "Band"
      trait_name <- colnames(dt)[t]
      predictor_cols <- grep("^Band_\\d", colnames(dt), value = TRUE)
      
      # Also handle possible extra space
      if (length(predictor_cols) == 0) {
        predictor_cols <- grep("^Band_\\s*\\d", colnames(dt), value = TRUE)
      }
      
      # Subset only complete rows for this trait + predictors
      dt_trait <- dt[, c(trait_name, predictor_cols)]
      dt_trait <- dt_trait[complete.cases(dt_trait), ]
      print(dim(dt_trait))
      print(colnames(dt_trait)[1:5])
      
      # --- Model 1: 5-fold CV model with selected ncomp ---
      cat("  Fitting CV model...\n")
      
      mod_cv <- plsr(
        as.formula(paste(trait_name, "~ .")), 
        data = dt_trait, 
        validation = "CV", 
        segments = 5
      )
      
      saveRDS(mod_cv, file = file.path(output_dir, paste0("mod_cv_", pt, "_", trait_name, "_iter",q,".rds")))
      
      # --- Model 2: Fit on KS calibration set with selected ncomp ---
      cat("  Fitting KS model...\n")
      mod_ks <- plsr(
        as.formula(paste(trait_name, "~ .")), 
        data = dt_trait[ken_mahal$model, ], 
        validation = "none"
      )
      
      saveRDS(mod_ks, file = file.path(output_dir, paste0("mod_ks_", pt, "_", trait_name, "_iter",q,".rds")))
    }
  }
}

# Run the simplified modeling step
FitAndSaveModels(pretreatment_results, ken_mahal_results, ncomp_overrides, traits = 1:5, output_dir = dir)

#{r Load Models Calc Performanc Metrics }
# Load and summarize saved PLSR model
SummarizeSavedModels <- function(pretreatment_results, ken_mahal_results, ncomp_overrides, traits = 1:5, input_dir = dir) {
  
  res_matr_list <- list()
  
  for (pt in names(pretreatment_results)) {
    cat("Summarizing Pretreatment:", pt, "\n")
    dt <- as.data.frame(pretreatment_results[[pt]])
    ken_mahal <- ken_mahal_results[[pt]]
    
    trait_names <- colnames(dt)[traits]
    res_matr <- matrix(NA, nrow = 8, ncol = length(trait_names))
    colnames(res_matr) <- trait_names 
    rownames(res_matr) <- c(
      "trait mean in all samples", "trait SD in all samples",
      "RMSEP in 5x CV in all samples", "r(pred,obs) in 5x CV in all samples",
      "RMSEP in KS test", "r(pred,obs) in KS test", 
      "NComp OneSigma","NComp used"
    )
    
    for (trait_name in trait_names) {
      task_id <- paste(pt, trait_name, sep = "_")
      
      cat("  Processing Trait:", trait_name, "\n")
      
      # Subset and filter
      t_index <- which(colnames(dt) == trait_name)
      dt_trait <- dt[, c(t_index, (length(trait_names)+1):ncol(dt))]
      dt_trait <- dt_trait[complete.cases(dt_trait), ]
      
      if (nrow(dt_trait) < 5) next
      
      # Load models
      mod_cv_path <- file.path(input_dir, paste0("mod_cv_", pt, "_", trait_name, "_iter",q,".rds"))
      mod_ks_path <- file.path(input_dir, paste0("mod_ks_", pt, "_", trait_name, "_iter",q,".rds"))
      
      if (!file.exists(mod_cv_path) || !file.exists(mod_ks_path)) {
        warning(paste("Missing model file(s) for", task_id))
        next
      }
      
      mod_cv <- readRDS(mod_cv_path)
      mod_ks <- readRDS(mod_ks_path)
      
      # Trait mean/SD
      res_matr["trait mean in all samples", trait_name] <- mean(dt_trait[, 1], na.rm = TRUE)
      res_matr["trait SD in all samples", trait_name] <- sd(dt_trait[, 1], na.rm = TRUE)
      
      # One-sigma component selection
      onesigma_ncomp <- as.numeric(selectNcomp(mod_cv, method = "onesigma", plot = FALSE))
      if (onesigma_ncomp < 1) onesigma_ncomp <- 1
      res_matr["NComp OneSigma", trait_name] <- onesigma_ncomp
      
      # User override
      ncomp.final <- ncomp_overrides[[pt]][trait_name]
      res_matr["NComp used", trait_name] <- ncomp.final
      
      # Safety: max comps in CV model
      max_comps_cv <- dim(mod_cv$validation$pred)[3]
      if (is.na(ncomp.final) || ncomp.final > max_comps_cv) {
        warning(sprintf("Reducing ncomp.final to %d for %s - %s (CV)", max_comps_cv, pt, trait_name))
        ncomp.final <- max_comps_cv
      }
      
      # CV metrics
      res_matr["RMSEP in 5x CV in all samples", trait_name] <- RMSEP(mod_cv, ncomp = ncomp.final)$val[[4]]
      res_matr["r(pred,obs) in 5x CV in all samples", trait_name] <- cor(mod_cv$validation$pred[,,ncomp.final],
                                                                         dt_trait[, 1], use = "pairwise.complete.obs")
      
      # KS prediction
      max_comps_ks <- ncol(mod_ks$scores)
      if (is.na(ncomp.final) || ncomp.final > max_comps_ks) {
        warning(sprintf("Reducing ncomp.final to %d for %s - %s (KS)", max_comps_ks, pt, trait_name))
        ncomp.final <- max_comps_ks
      }
      
      pred_test <- predict(mod_ks, ncomp = ncomp.final, newdata = dt_trait[-ken_mahal$model, ])
      res_matr["RMSEP in KS test", trait_name] <- sqrt(mean((dt_trait[-ken_mahal$model, 1] - pred_test)^2, na.rm = TRUE))
      res_matr["r(pred,obs) in KS test", trait_name] <- cor(as.vector(pred_test), dt_trait[-ken_mahal$model, 1], use = "pairwise.complete.obs")
      
      # Optional cleanup
      #rm(mod_cv, mod_ks); gc() #CHD commented for diagnostic purposes
    }
    
    res_matr_list[[pt]] <- res_matr
  }
  
  return(res_matr_list)
}

# Call the function
final_results <- SummarizeSavedModels(pretreatment_results, ken_mahal_results, ncomp_overrides, traits = 1:5, input_dir = dir)#change to 1:8 when adding starch

# Combine all data frames in the list `final_results` with unique row names
final_combined <- do.call(rbind, lapply(names(final_results), function(name) {
  df <- final_results[[name]]
  # Add the prefix to the row names
  rownames(df) <- paste0(name, "_", rownames(df))
  return(df)
}))

write.csv(final_combined, paste0("NIRS_custom_calibrations_wPretreats_PLSR_barley_results_", Sys.Date(),"_iter",q,".csv")) 
} #end for 1 in niter

#Before proceeding review results file output by previous chunk******
# This chunk loads in fitted models selects best and predicts over the full set
setwd("saved_models/2026-01-09") 
mod <- readRDS("mod_cv_snv.sg2_Ash.rds")

# List model files
mods <- list.files()

#Will use iter5
#Choose best model using the output from the last chunk************
best_mods <- mods[grepl("^mod_cv_snv\\.sg_.*\\_iter5.rds$", mods)]

# Load models into a named list
loaded_models <- list()
for (modfile in best_mods) {
  model_name <- sub("\\_iter5.rds$", "", modfile)  # remove .rds extension
  loaded_models[[model_name]] <- readRDS(modfile)
}

# Create a named list to store predictions
predictions <- list()

# Loop through each loaded model and predict
for (i in names(loaded_models)) {
  tryCatch({
    pred_full <- predict(loaded_models[[i]], snv.sg)  # **** Choose best pretreatment
    
    # Determine trait name
    trait_name <- sub("mod_cv_snv\\.sg_", "", i)
    
    # Find how many components were selected
    ncomp_selected <- ncomp_overrides[["snv.sg"]][trait_name]
    
    print(paste("Using", ncomp_selected, "components to predict", trait_name, "." ))
    
    if (is.null(ncomp_selected)) {
      warning(glue("No ncomp override found for {trait_name}, defaulting to 1 component."))
      ncomp_selected <- 1
    }
    
    # Slice the prediction to the selected number of components
    predictions[[trait_name]] <- pred_full[, , ncomp_selected]
    
  }, error = function(e) {
    warning(glue("Prediction failed for {i}: {e$message}"))
    predictions[[trait_name]] <- NA
  })
}

# Turn the list into a dataframe
pred_df <- as.data.frame(predictions)

#Make column for merging
pred_df$sample_number <- make.unique(names(predictions[[1]]))

write.csv(pred_df,"pred_df_from20260109models_iter5.csv",quote=F,row.names = F)
#Set output dir
#output_dir <- "Custom Calibration/plots" 
#today <- Sys.Date()

#dat_wetchem$Moisture <- 100 - dat_wetchem$DM

#remove Ntotal column from all DF
#wetchem <- dat_wetchem
#wetchem= wetchem[, !colnames(wetchem) %in% "Ntotal"]

whole_macros_cleaned = BarleyF2.whole.macs[, c(5,9,13,17, 25)]
colnames(whole_macros_cleaned) <- gsub("\\..[2]+\\.", "",colnames(whole_macros_cleaned))

#pred_df_cln= pred_df[, !colnames(pred_df) %in% "Ntotal_Pred"]

# Standardize pred_df + wetchem + whole grain values column names
if (names(pred_df)[1] != "Ash_Pred") {
  names(pred_df)[1:5] <- paste0(names(pred_df)[1:5], "_Pred")
}
if (names(dat_wetchem_NoOuts)[2] != "Ash_WC") {
  names(dat_wetchem_NoOuts)[2:6] <- paste0(names(dat_wetchem_NoOuts)[2:6], "_WC")
}
if (names(BarleyF2.whole.macs)[5] != "Ash_CG") {
  names(whole_macros_cleaned)[2:5] <- paste0(names(whole_macros_cleaned)[2:5], "_CG")#change to 1:8 when adding starch
}

#change wetchem "plot" and whole_grain_macros "sample.number" column to sample_number
names(dat_wetchem_NoOuts)[names(dat_wetchem_NoOuts) == 'plot'] <- 'sample_number'
names(whole_macros_cleaned)[names(whole_macros_cleaned) == 'Sample.Number'] <- 'sample_number'

# Merge
CG_vs_Cust_vs_WetChem <- dat_wetchem_NoOuts %>%
  inner_join(pred_df, by = "sample_number") %>%
  inner_join(whole_macros_cleaned, by = "sample_number")

# Trait map
traits <- list(
  Ash = list(pred_col = "Ash_Pred", CG_col = "Ash_CG", wetchem_col = "Ash_WC"),
  Fat = list(pred_col = "Fat_Pred",CG_col = "Fat_CG", wetchem_col = "Fat_WC"),
  Moisture = list(pred_col = "Moisture_Pred", CG_col = "Moisture_CG", wetchem_col = "Moisture_WC"),
  Protein = list(pred_col = "Protein_Pred", CG_col = "Protein_CG", wetchem_col = "Protein_WC"))

# R² calculator
calc_r2 <- function(x, y) {
  if (is.null(x) || is.null(y)) return(NA)
  df <- na.omit(data.frame(x = x, y = y))
  if (nrow(df) == 0) return(NA)
  round(cor(df$x, df$y)^2, 3)
}

# Global limits
wetchem_cols <- sapply(traits, function(x) x$wetchem_col)
pred_cols    <- sapply(traits, function(x) x$pred_col)
CG_cols     <- sapply(traits, function(x) x$CG_col)

existing_wetchem_cols <- intersect(wetchem_cols, names(CG_vs_Cust_vs_WetChem))
global_xrange <- range(unlist(CG_vs_Cust_vs_WetChem[existing_wetchem_cols]), na.rm = TRUE)

existing_pred_cols <- intersect(pred_cols, names(CG_vs_Cust_vs_WetChem))
existing_CG_cols  <- intersect(CG_cols,  names(CG_vs_Cust_vs_WetChem))

global_yrange <- range(unlist(CG_vs_Cust_vs_WetChem[c(existing_pred_cols, existing_CG_cols)]), na.rm = TRUE)

for (trait in names(traits)) {
  pred_col <- traits[[trait]]$pred_col
  CG_col  <- traits[[trait]]$CG_col
  wet_col  <- traits[[trait]]$wetchem_col
  x <- CG_vs_Cust_vs_WetChem[[wet_col]]
  
  # --- Custom vs WetChem ---
  y_custom <- CG_vs_Cust_vs_WetChem[[pred_col]]
  r2_custom <- calc_r2(x, y_custom)
  
  if (!all(is.na(x)) && !all(is.na(y_custom))) {
    p1 <- ggplot(CG_vs_Cust_vs_WetChem, aes(x = .data[[wet_col]], y = .data[[pred_col]])) +
      geom_point(alpha = 0.7, color = "blue", size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      annotate("text", x = min(x, na.rm = TRUE), y = max(y_custom, na.rm = TRUE),
               label = paste0("R² = ", r2_custom), hjust = 0, vjust = 1.2, size = 8) +
      labs(
        title = paste(trait, "- Custom Calibration vs Wet Chemistry"),
        subtitle = "Custom calibration: PLSR with SNV and SG",
        x = "Wet Chemistry (%)",
        y = "Custom Prediction (%)"
      ) +
      theme_bw(base_size = 22) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = margin(1, 1, 1, 1, "cm"))+
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      scale_y_continuous(labels = function(y) sprintf("%.1f", y))
    
    ggsave(
      filename = file.path(paste0("Custom_SNV.SG_vs_WetChem_", trait, "_", Sys.Date(), ".png")),
      plot = p1,
      width = 10, height = 10  # Consistent PDF size
    )
    
    # --- CG vs WetChem ---
    if (!is.na(CG_col)) {
      y_CG <- CG_vs_Cust_vs_WetChem[[CG_col]]
      r2_CG <- calc_r2(x, y_CG)
      
      p2 <- ggplot(CG_vs_Cust_vs_WetChem, aes(x = .data[[wet_col]], y = .data[[CG_col]])) +
        geom_point(alpha = 0.7, color = "darkgreen", size = 3) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        annotate("text", x = min(x, na.rm = TRUE), y = max(y_CG, na.rm = TRUE),
                 label = paste0("R² = ", r2_CG), hjust = 0, vjust = 1.2, size = 8) +
        labs(
          title = paste(trait, "- Manufacturer Calibration vs Wet Chemistry"),
          subtitle = "Manufacturer calibration: PLSR or ANN",
          x = "Wet Chemistry (%)",
          y = "Manufacturer Prediction (%)"
        )+
        theme_bw(base_size = 22) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.margin = margin(1, 1, 1, 1, "cm"))+
        scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
        scale_y_continuous(labels = function(y) sprintf("%.1f", y))
      
      ggsave(
        filename = file.path(paste0("CG_vs_WetChem_", trait, "_", Sys.Date(), ".png")),
        plot = p2,
        width = 10, height = 10
      )
    }
  }
}

#{r write final results file}

# Merge predictions and CG columns
CG_vs_Cust <- merge(pred_df, whole_macros_cleaned, by = "sample_number")

# Merge wet chemistry into final output
CG_vs_Cust_plusWet <- merge(CG_vs_Cust, dat_wetchem, by = "sample_number", all.x = TRUE)

# Write final results
write.csv(
  CG_vs_Cust_plusWet,
  paste0("CustPreds_Barley_F2_SNV.SG_", Sys.Date(), ".csv"),
  row.names = FALSE
)
