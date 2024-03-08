
##### Mapping the bulk expression matrix #####

BiocManager::install("umap")
BiocManager::install("uwot")

library(uwot)
library(ggplot2)

# Load the tirosh data from *DTD tutorials.R* file
tirosh_data <- readRDS("Data/DTD tutorials.rds")
bulk.expres <- tirosh_data$bulk_expres

##### Dimentionality reduction using uwot package #####

# We are going to use Uniform Manifold Approximation function
# This UMAP function allows more complicated distributions
# that it learns from the data.
# The advantage of this feature is that UMAP can do a better job 
# separating clusters, especially when some of those clusters 
# may be more similar to each other than others.

View(bulk.expres)


set.seed(142)   #To get reproducible expression value

umap.data <- umap(t(bulk.expres),n_neighbors = 3, learning_rate = 0.5, 
                  init = "random", n_epochs = 1000, metric = "correlation")

# Analyse the umap data. default n_components (dimension) is 2.
View(umap.data)

#Convert matrix into data frame to plot
data.plot <- data.frame(UMAP1  = umap.data[,1],
                        UMAP2 = umap.data[,2], 
                        tumors = colnames(bulk.expres))
                  
View(data.plot)

#Plotting the umap data
my.plot <- ggplot(data = data.plot,
                   mapping = aes( x = UMAP1, y = UMAP2 )) 
my.plot +
 geom_point(color = colnames(bulk.expres))+
  labs( x = "UMAP1", y = "UMAP2", subtitle = "bulk.expres UMAP",  )

ggsave("UMAP_plot_bulk_expres.tiff",path = "Results", units="in", width=5, height=4, dpi=300, compression = 'lzw')



