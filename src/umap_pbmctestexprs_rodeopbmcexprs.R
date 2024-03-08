
#### Create expression matrix using pbmc reference profile and cell proportion matrix ####

# Load the data from *cor_estcwithDTD_pbmcpheno.R* file
pbmc_data2 <- readRDS("Data/rodeo_pbmc_data")
ref_profile_pbmc <- pbmc_data2$ref_rodeo_pbmc
test_pheno_pbmc <- pbmc_data2$pbmc_test_pheno
test_bulk_pbmc <- pbmc_data2$pbmc_test_expres


# Y = X.C (Where, X is reference profile, C is cell-proportion matrix and Y = bulk expression matrix)

# Create bulk expression matrix by multiplying X and C matrices

estimated_bulk_pbmc <- ref_profile_pbmc %*% test_pheno_pbmc

# UMAP with rodeo generated bulk expression profile and test bulk expression 

merge_mat_pbmc<- cbind(test_bulk_pbmc, estimated_bulk_pbmc)

# Reduce dimentionality of combined matrix using UMAP function

library(uwot)
library(ggplot2)

set.seed(123)

umap_data_pbmc <- umap(t(merge_mat_pbmc),n_neighbors = 7, n_components = 3, learning_rate = 1, 
                  init = "random", n_epochs = 1000, metric = "correlation")

# Create vectors to make umap dataframe
group_shape <- c(rep("Test bulk", ncol(test_bulk_pbmc)), rep("Estimated bulk",ncol(estimated_bulk_pbmc)))
group_color <- colnames(merge_mat_pbmc)

# Create a dataframe
umap_plot_data_pbmc <- data.frame(UMAP1 = umap_data_pbmc[,1],
                                  UMAP2 = umap_data_pbmc[,2],
                                  group = group_shape,
                                  color = group_color)



# Use ggplot to plot the UMAP
ggplot(data = umap_plot_data_pbmc, mapping = aes(UMAP1, UMAP2))+
  geom_point(aes(shape = as.factor(group), color = group_color))+
  scale_shape_manual(name= "Bulk expression matrix", values = c(17, 15) ,guide = "legend") +
  scale_color_manual( values= rainbow(40), guide = "none")+
  theme(legend.position = "right")+
  labs(subtitle = "PBMC bulk expression data")
  
ggsave("UMAP_plot731c_pbmcdata.tiff", path = "Results", units="in", width=5, height=4, dpi=300, compression = 'lzw')
