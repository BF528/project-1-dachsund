# install packages
install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("sva")
BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2.db")

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)

###### set working directory to where the .CEL files are 
setwd("/projectnb/bf528/users/dachshund/project_1/samples/all_samples")
# getwd()


###### robust multi-array average (RMA) expression ----

# read .CEL files in current directory and normalize them
# cel_path = '/projectnb/bf528/users/dachshund/project_1/samples'
data = ReadAffy()
norm_data <- rma(dat)

# get RMA normalized expression values for each array (133 samples, 54675 features)
# the ExpressionSet class combines different pieces of into into a single convenient structure
# data includes expression data from microarrays , gene annotations
show(norm_data)
rma_values <- exprs(norm_data)
rma_values_df <- as.data.frame(rma_values)


###### fitPLM ----

# Fit A Probe Level Model To Affymetrix Genechip Data.
# fitPLM converts an AffyBatch into an PLMset by fitting a specified robust linear model to the probe level data.

pset <- fitPLM(data, normalize=TRUE, background=TRUE)
(pset)


##### RLE ----
# compute the relative log expression (RLE) of the microarray samples

# RLE for each sample
RLE_scores <- RLE(pset, type="stats")

# plot median RLE in a boxplot
RLE(pset, main="RLE")

# plot median RLE in a histogram
hist(RLE_scores["median",], xlab = "Median RLE", ylim = c(0, 60), col = c("darkred"), main = NULL)


# NUSE ----
# the standard error estimates obtained for each gene on each array from fitPLM
# are taken and standardized across arrays so that the median standard error for that
# genes is 1 across all arrays. This process accounts for differences in variability between
# genes. An array where there is elevated SE relative to the other arrays is typically of
# lower quality. 

# compute normalized unscaled standard error (NUSE)
NUSE_scores <- NUSE(pset, type= "stats")

# plot NUSE in a histogram
hist(NUSE_scores["median",], xlab = "Median NUSE", ylim = c(0, 50), col = c("lightblue"), main = NULL)


# ComBat ----
# use ComBat to correct for batch effects. ComBat allows users to adjust for 
# batch effects in datasets where the batch covariate is known, using methodology
# described in Johnson et al. 2007. It uses either parametric or non-parametric
# empirical Bayes frameworks for adjusting data for batch effects. 
# Users are returned an expression matrix that has been corrected for batch effects. 

# transform normalized data into a matrix in which rows are probesets and columns are samples
norm_data_matrix <- exprs(norm_data)
dim(norm_data_matrix)

# path to annotation file on SCC, contains a clinical and batching annotation used 
# by the authors for their analysis
meta_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")

# Batch effects include both Center and RNA extraction method and have been merged
# into a single variable called *normalizationcombatbatch* in the annotation file
bat = meta_data$normalizationcombatbatch
dim(bat)

# mod is a model matrix that contains biological covariates, including the outcome of 
# interest. This ensures that biological variability is preserved. Our features
# of interest are tumor and MMR status. They have been merged into a single 
# variable called *normalizationcombatmod*
model <- model.matrix(~normalizationcombatmod, data = meta_data)
dim(model)

# correct for batch effects while preserving features of interest - tumor and MMR status
combat_corrected <- ComBat(dat = norm_data_matrix, batch = bat, mod = model)


###### principal component analysis (PCA) ----
# perform PCA on the batch corrected and normalized data
# PCA finds a new coordinate system for multivariate data such that the first coordinate has maximal variance, 
# the second coordinate has maximal variance subject to being orthogonal to the first

# center and scale data (this is done within a column and within each gene rather than sample, hence the need to transpose)
transposed_combat_corrected <- t(combat_corrected)
scaled_transposed <- scale(transposed_combat_corrected)

# return the scaled transposed data back to its original orientation
untransposed_scaled_data <- t(scaled_transposed)

# perform a principal components analysis on the given data matrix 
# the prcomp function provides a resulting object that gives us standard deviation,
# proportion of variance explained by each principal component, 
# and the cumulative proportion of variance explained
# scale = false means that data is not scaled
# center = false means that input variables are not shifted to be zero centered
pca <- prcomp(untransposed_scaled_data, scale = FALSE, center = FALSE)
colnames(pca)

# view the values for each of the principal components by accessing the $rotation attribute of your prcomp object.
pca_rot <- pca$rotation
pca_rot_df <- as.data.frame(pca_rot)
head(pca$rotation)
pc_values <- pca$rotation[,1:2]

# determine the principal components
PC1 <- pca$rotation[,1] # the first principal component
PC2 <- pca$rotation[,2] # the second principal component 

# examine the percent variability explained by each principal component by looking at the $importance attribute
pca_summary <- summary(pca)
(pca_summary$importance)

# plot PC1 vs PC2 with the percent variability attributed to these principal components shown
plot(PC1, PC2, 
     xlab = paste0("PC1 (", round(pca_summary$importance[2,1]*100, digits = 2), "%)"), 
     ylab = paste0("PC2 (", round(pca_summary$importance[2,2]*100, digits = 2), "%)"),
     cex = 1,
     main = "PC1 vs PC2")

# plot PC1 vs PC2 with the percent variability using ggplot (just showing another method of plotting)
pc_df <- as.data.frame(pc_values)
ggplot(pc_df, mapping = aes(x = PC1, y = PC2)) + 
  geom_point(size=1.5) +
  # labs(title = "PCA plot") +
  xlab(paste0("PC1:",format(pca_summary$importance[2,1]*100, digits = 4), "%")) +
  ylab(paste0("PC2:", format(pca_summary$importance[2,2]*100, digits = 3), "%"))

# looking for outliers in PC1 and PC2 are greater than 3 standard deviations
pc_outliers <- which(PC1 > mean(PC1) + 3*sd(PC1) | PC1 < mean(PC1) - 3*sd(PC1) | 
                       PC2 > mean(PC2) + 3*sd(PC2) | PC2 < mean(PC2) - 3 * sd(PC2)) 
(pc_outliers)

# another way of showing outliers using the previously stated parameters
pc_outliers_id <- which(pca_rot_df$PC1 > mean(pca_rot_df$PC1) + 3*sd(pca_rot_df$PC1) | pca_rot_df$PC1 < mean(pca_rot_df$PC1) - 3*sd(pca_rot_df$PC1) |
                          pca_rot_df$PC2 > mean(pca_rot_df$PC2) + 3*sd(pca_rot_df$PC2)| pca_rot_df$PC2 < mean(pca_rot_df$PC2) - 3*sd(pca_rot_df$PC2))
(pc_outliers_id)

# plot outliers in PC1 and PC2
ggplot(data = pc_df, mapping = aes(y = PC1)) + 
  geom_boxplot() + 
  annotate('text', label='A', x=-Inf, y=Inf, hjust=-0.5, vjust=1.5, size = 6) + 
  labs(x = " ")

ggplot(data = pc_df, mapping = aes(y = PC2)) + 
  geom_boxplot(outlier.size = 2) + 
  annotate('text', label='B', x=-Inf, y=Inf, hjust=-0.5, vjust=1.5, size = 6) + 
  labs(x = " ")

# check outliers for PC2
which(PC2 > mean(PC2) + 3*sd(PC2) | PC2 < mean(PC2) - 3*sd(PC2)) 

# remove outliers from the ComBat adjusted data
id <- which(!(pca_rot_df$PC1 > mean(pca_rot_df$PC1) + 3*sd(pca_rot_df$PC1) | pca_rot_df$PC1 < mean(pca_rot_df$PC1) - 3*sd(pca_rot_df$PC1) |
                pca_rot_df$PC2 > mean(pca_rot_df$PC2) + 3*sd(pca_rot_df$PC2)| pca_rot_df$PC2 < mean(pca_rot_df$PC2) - 3*sd(pca_rot_df$PC2)))
(id)
dim(id)

# ComBat adjusted data with outlier(s) removed 
combat_corrected_no_outliers <- combat_corrected[, id]
head(combat_corrected_no_outliers)
dim(combat_corrected_no_outliers)

# plot PC1 vs PC2 with outliers removed
pca_no_outliers <- prcomp(untransposed_scaled_data[,id], scale = FALSE, center = FALSE)
pc_values_no_outliers <- pca_no_outliers$rotation[,1:2]
pc_no_outliers_df <- as.data.frame(pc_values_no_outliers)
pca_no_outliers_summary <- summary(pca_no_outliers)
pca_no_outliers_summary$importance
ggplot(pc_no_outliers_df, mapping = aes(x = PC1, y = PC2)) + 
  geom_point(size=1.5) +
  # labs(title = "PCA plot") +
  annotate('text', label='C', x=-Inf, y=Inf, hjust=-0.5, vjust=1.5, size = 6) +
  xlab(paste0("PC1:",format(pca_no_outliers_summary$importance[2,1]*100, digits = 4), "%")) +
  ylab(paste0("PC2:", format(pca_no_outliers_summary$importance[2,2]*100, digits = 3), "%")) 



##### write csv files ----
setwd("/projectnb/bf528/users/dachshund/project_1")

# write RMA-normalized gene expression values to a file
write.csv(rma_values, file = "rma.csv")
rma_file <- read.csv(file = 'rma.csv')
head(rma_file)

# write ComBat expression data (that is also RMA normalized) to a CSV file in which probesets are rows and samples are columns
write.csv(combat_corrected, file = "ComBat_adjusted_data.csv")

# write ComBat adjusted data without outlier(s)
write.csv(combat_corrected_no_outliers, file = "ComBat_adjusted_outliers_removed_data.csv")


