######################################################
#                        k-means                     #
######################################################

# simulate normally distributed data and store it in a matrix (100 samples and 2 variables)
x <- matrix(rnorm(n=100*2), nrow=100, ncol=2)

## Artificially create 4 clusters
# simulate the 4 cluster means
xmean <- matrix(rnorm(n=8, sd=5), nrow=4, ncol=2)
# assign each sample to one of the four clusters
ind <- sample(1:4, 100, replace=TRUE)
# add the means of their corresponding cluster
x <- x + xmean[ind,]

# Plot the dataset we created 
# the data point color shows the clusters we simulated
plot(x, col=ind, pch=19)
# the data point shape shows the clusters we simulated
plot(x, pch=ind)

# We know the "true" cluster IDs, but we won't give it as an input to the kmeans algorithm.
# Check the help of the kmeans function
?kmeans

# Run kmeans asking for k=3 clusters and 20 random initializations. 
km3 <- kmeans(x, centers=3, nstart=20)
# Run kmeans asking for k=4 clusters and 20 random initializations. 
km4 <- kmeans(x, centers=4, nstart=20)
# Run kmeans asking for k=4 clusters and only one 1 random initialization. 
km41 <- kmeans(x, centers=4, nstart=1)

# Look at the kmeans output
km4
# cluster assignment for each individual
km4$cluster
# within-cluster sum of squares
km4$withinss
# The number of individuals in each cluster
km4$size
# The means of the clusters (first and second variable)
km4$centers

# Plot the kmean output for k=3
plot(x, col=km3$cluster, pch=19)
# Plot the kmean output for k=4 and 20 random initializations
plot(x, col=km4$cluster, pch=19)
# Plot the kmean output for k=4 and 1 random initialization
plot(x, col=km41$cluster, pch=19)

######################################################
#                      k-medoids                     #
######################################################

# Load library cluster
library(cluster)
# Run PAM algorithm on the data we simulated asking for 3 or 4 clusters
pam3 <- pam(x, 3)
pam4 <- pam(x, 4)
#  Plot the PAM output for four clusters
plot(x, col=pam4$cluster, pch=19)
#  Circle the data points that are the medoids of each cluster
points(pam4$medoids, cex=2, col="purple", lwd=2)

par(mfrow = c(1,2))    # to plot the two plots next to each other
plot(silhouette(pam3)) # plot silhouette of PAM with k=3
plot(silhouette(pam4)) # plot silhouette of PAM with k=4

library(factoextra)
fviz_nbclust(x, pam, method="silhouette")
fviz_nbclust(x, pam, method="gap_stat")

# Load the example dataset
load("dat_pca_clus_counts.rda")
# use log2 expression
dt[,3:ncol(dt)] <- log2(dt[,3:ncol(dt)])
# Remove the first two annotation columns
dt <- dt[, c(-1,-2)]
# dt is ready for the exercise!

library(cluster)
res.pam <- pam(dt, k=3)
res.pam2$clustering

library(factoextra)
fviz_nbclust(dt, FUNcluster=pam, method="silhouette")
# generates an error as trying to extract more clusters than individuals
fviz_nbclust(dt, FUNcluster=pam, method="silhouette", k.max=5)

res.pam2 <- pam(dt, k=2)
res.pam2$clustering

######################################################
#               hierarchical clustering              #
######################################################
# Load again the example dataset, this time we keep it on the read count scale
load('dat_pca_clus_counts.rda')

# Extract and remove the two annotation columns in order to only keep the quantitative variables
annot <- dt[,c(1,2)]
dt <- dt[,-c(1,2)]

# Run hierarchical clustering using the hclust() function with different linkages and different transformations of the data.
# the function dist() computes the pairwise euclidean distance
# on the count scale
eucl = dist(dt)
hc.complete <- hclust(eucl, method = 'complete')
hc.average <- hclust(eucl, method = 'average')
hc.single <- hclust(eucl, method = 'single')
hc.ward <- hclust(eucl, method = 'ward.D')

# on the log2(count) scale
eucllog <- dist(log2(dt))
hc.completelog <- hclust(eucllog, method = 'complete')
hc.averagelog <- hclust(eucllog, method = 'average')
hc.singlelog <- hclust(eucllog, method = 'single')
hc.wardlog <- hclust(eucllog, method = 'ward.D')

# Plot the dendrograms resulting from each method
par(mfrow=c(2,2))
plot(hc.average, main = 'Average linkage \n no log transformation')
plot(hc.averagelog, main = 'Average linkage \n log transformation')
plot(hc.ward, main = 'Ward  \n no log transformation')
plot(hc.wardlog, main = 'Ward  \n log transformation')

######################################################
#                       heatmaps                     #
######################################################
# We will plot a heatmap of our dataset with a dendrogram for both genes and samples

# load packages
library(RColorBrewer) ## package useful for colors
library(gplots)       ## package containing a heatmap function among other things

# Choose a color for each patient
colpatients = c(p1='red', p2='turquoise', p3='blue', p4='green')

# Choose a range of colors for the heatmap
# display the color palettes available in the package
display.brewer.all()
# choose for example the "spectral" palette and create a gradient of colors
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# We will use euclidean distance and Ward method to build the dendrograms
hclusfun <- function(x) hclust(dist(x), method="ward.D")

## Heatmap without centering or scaling the genes 
heatmap.2(t(dt), dendrogram="both", scale="none",
          ColSideColors=colpatients[annot$patient], 
          trace="none", key=TRUE, col=myPalette, hclustfun = hclusfun)
## Heatmap while centering the genes 
heatmap.2(t(dt), dendrogram="both", scale="row", 
          ColSideColors=colpatients[annot$patient], trace="none", col=myPalette, 
          hclustfun = hclusfun)

######################################################
#                       exercise                     #
######################################################
data(iris)

# hierarchical clustering with euclidean distance (dist's default)
hc_eucl_av <- hclust(d = dist(iris[,1:4]), method = 'average')
hc_eucl_compl <- hclust(d = dist(iris[,1:4]), method = 'complete')
# hierarchical clustering with 1 - Pearson correlation as distance
dcor <- 1 - cor(t(iris[,1:4]))
hc_cor_av <- hclust(d = as.dist(dcor), method = 'average')
hc_cor_compl <- hclust(d = as.dist(dcor), method = 'complete')

# plot the dendrograms with species as labels 
par(mfrow = c(2,2))
plot(hc_eucl_av, labels = iris[,5], main = 'Euclidean distance, average linkage')
plot(hc_eucl_compl, labels = iris[,5], main = 'Euclidean distance, complete linkage')
plot(hc_cor_av, labels = iris[,5], main = '1-cor distance, average linkage')
plot(hc_cor_compl, labels = iris[,5], main = '1-cor distance, complete linkage')

# function that will be applied to each clustering result
# cut the dendograms into 3 clusters 
# and print the number of species representative in each
cutree_hc <- function(res.hc, species, label) {
  clusters <- cutree(res.hc, 3)
  table(clusters, species)
}

# function that will be applied to each clustering result
cutree_hc(hc_eucl_av, species = iris[,5])
cutree_hc(hc_eucl_compl, species = iris[,5])
cutree_hc(hc_cor_av, species = iris[,5])
cutree_hc(hc_cor_compl, species = iris[,5])

# heatmap
library(gplots)
library(RColorBrewer)
hclusfun <- function(x) hclust(x, method = 'average') 
distfun <- function(x) as.dist(1-cor(t(x)))
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
heatmap.2(as.matrix(iris[,1:4]), scale="none", 
          RowSideColors=c(setosa = 'blue', versicolor = 'pink', virginica = 'red')[as.character(iris[,5])],
          trace="none", key=TRUE, col=myPalette, 
          hclustfun = hclusfun, distfun = distfun)
