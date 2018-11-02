## load the factoextra and FactoMineR packages
## (install first if necessary)
library(factoextra)
library(FactoMineR)

## Set your working directory location
setwd("~/Teaching/Cours_PhD/PCA/")

## Load example dataset
load('dat_pca_clus_counts.rda')


## One data.frame object has been loaded:
## RNASeq read count per gene, data restricted to 50 genes
head(dt)
dim(dt)

## PCA will be computed on the log2 read counts
dt[,3:ncol(dt)] <- log2(dt[,3:ncol(dt)])

## Run the PCA 
respca = PCA(dt, scale.unit=TRUE, graph=F, ncp=5, quali.sup=1:2)

# X: the data set used. Rows are individuals and columns are 	numeric variables.
# scale.unit: if TRUE, the data are scaled to unit variance 			 before the analysis. 
# ncp: number of PCs to keep.
# graph: if TRUE, the graphs are displayed.

respca

# The attribute eig contains information relative to the variance explained in each component.
eig = get_eig(respca)

# Plot the % of variances explained by each PC
fviz_screeplot(respca, addlabels = TRUE)

##  Extract the data associated with the observations from the pca object
ind = get_pca_ind(respca)
ind

## coordinates of the individuals in the PCA space
ind$coord


## Plot of the samples onto the space spanned by the first and the second PCs.
## Color the data points according to the patient

fviz_pca_ind(respca, 
             axes = c(1,2), 
             habillage = 'patient', 
             invisible = 'quali')

# Equivalent with the function plot.PCA (run faster for bigger datasets)
plot.PCA(respca, axes = c(1,2), choix = 'ind', habillage = 1,  invisible = 'quali')

## Plot of the samples onto the space spanned by the third and the fourth PCs.
## Color the data points according to the patient

fviz_pca_ind(respca, 
             axes = c(3,4), 
             habillage = 'patient', 
             invisible = 'quali')

# Equivalent with the function plot.PCA (run faster for bigger datasets)
plot.PCA(respca, axes = c(3,4), choix = 'ind', habillage = 1,  invisible = 'quali')

## Plot of the samples onto the space spanned by the third and the fourth PCs.
## Color the data points according to the treatment

fviz_pca_ind(respca, 
             axes = c(3,4), 
             habillage = 'treatment', 
             invisible = 'quali')

# Equivalent with the function plot.PCA (run faster for bigger datasets)
plot.PCA(respca, axes = c(3,4), choix = 'ind', habillage = 2,  invisible = 'quali')

##  Extract the data associated with the observations from the pca object
ind = get_pca_ind(respca)
ind

## squared cosines of the observations
ind$cos2

##  Extract the data associated with the variables from the pca object
var = get_pca_var(respca)
var

ind$contrib
colSums(ind$contrib)

fviz_pca_var(respca, col.var = "cos2", gradient.cols = c("orange", "blue"))

## Plot of the correlation of the variables with PC3 and PC4.

fviz_pca_var(respca, axes = c(3,4), choix = 'var')

## You can visualize the cos2 of the variables on all the dimensions using the corrplot package

library("corrplot")
corrplot(respca$var$cos2)

ibrary(factoextra) 
fviz_pca_biplot(respca, axes = c(1,2), 
                habillage = 'patient', 
                select.var = list(cos2 = 0.99),
                invisible = c("ind.sup", "quali"))

## Plot of the correlation of the variables with the PCs
## use parameter lim.cos2.var to only show variables with a cos2 value larger than 0.5

fviz_pca_var(respca, axes = c(3,4), choix = 'var', select.var = list(cos2 = 0.5))



## Plot of the correlation of the variables with the PCs
## use parameter lim.cos2.var to only show variables with a cos2 value larger than 0.5

fviz_pca_var(respca, axes = c(3,4), choix = 'var', select.var = list(cos2 = 0.5))

## You can visualize the cos2 of the variables on all the dimensions using the corrplot package

library("corrplot")
corrplot(respca$var$cos2)

library(factoextra) 
fviz_pca_biplot(respca, axes = c(1,2), 
                habillage = 'patient', 
                select.var = list(cos2 = 0.99),
                invisible = c("ind.sup", "quali"))

fviz_pca_biplot(respca, axes = c(3,4), 
                habillage = 'treatment', 
                select.var = list(cos2 = 0.5),
                invisible = c("quali"))

## Declare the two samples from patient 3 as supplementary individuals.
respca = PCA(X = dt, scale.unit = TRUE, graph = F, 
             ncp = 5, quali.sup = 1:2,
             ind.sup = which(dt$patient == 'p3'))
fviz_pca_ind(respca, axes = c(1,2), 
             habillage = 'patient’,
	    invisible = c("ind.sup", "quali"))


## Project the patient 3 samples in the PC space defined using only patients 1 and 2
fviz_pca_ind(respca, axes = c(1,2), 
                        habillage = 'patient’,
             invisible = c("quali"))

