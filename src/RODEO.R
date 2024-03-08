##### RODEO #####

devtools::install_github("elolab/Rodeo", force = TRUE)

# Load required packages
library(MASS)

library(Rodeo)

# Load training data from *DTD tutorials.R* file

tiroshdata <- readRDS(file = "Data/DTD tutorials.rds")
training.data <- tiroshdata$train
testing.data <- tiroshdata$test
X.matrix <- tiroshdata$X_matrix
bulk.expres <- tiroshdata$bulk_expres


# Rodeo package needs bulk expression and cellular composition as input.

reference_x <- Rodeo(training.data$mixtures, training.data$quantities)


#The input are bulk expression and cell composition from training data.
#This provides reference matrix with genes in rows and cell types in column.
#reference matrix provides how strongly each cell type expressed in each gene.

library(DTD)
set.seed(1234)
library(ggplot2)

# Estimating cell-composition using test-bulk & estimated reference matrix

estimate_c_ref <- estimate_c(X.matrix = reference_x,
                              new.data = testing.data$mixtures,
                              DTD.model = rep(1,nrow(reference_x)),
                              estimate.c.type = "direct")

# Correlation between test and estimated_ref cell type composition matrices

cor_cell_ref_test <- mapply(cor, as.data.frame(t(estimate_c_ref)), 
                            as.data.frame(t(testing.data$quantities)))
cor_cell_ref_test

# This provides correlation between each cell type 
# Across all mixture-samples from two matrices.


# Estimating cell-composition from test bulk and X.matrix

estimate_c_X <- estimate_c(   X.matrix = X.matrix,
                              new.data = testing.data$mixtures,
                              DTD.model = rep(1,nrow(X.matrix)),
                              estimate.c.type = "direct")

# Correlation between test and estimated_X cell type composition matrices

cor_cell_X_test <- mapply(cor, as.data.frame(t(estimate_c_X)), 
                            as.data.frame(t(testing.data$quantities)))
cor_cell_X_test




####  RODEO Estimated and generated reference profile analysis ####

create_plots <- function(cell_type, reference_x, X.matrix) {
  
  data_cell <- data.frame(RODEO = reference_x[, cell_type],
                          Gen.Ref = X.matrix[, cell_type],
                          Celltype = cell_type)
  
  cell_plot <- ggplot(data = data_cell,
                      mapping = aes(x = RODEO, y = Gen.Ref)) +
    geom_point(aes(color = Celltype)) +
    scale_color_manual(values = c("B" = "blue", "CAF" = "red", "Macro" = "black", "NK"= "orange", "T" = "purple")) +
    labs(subtitle = paste0(cell_type, "_plot"),
         x = "RODEO generated reference matrix ",
         y = "X.matrix reference matrix")
  
  ggsave(paste0(cell_type, "_plot.tiff"),
         path = "/Results",
         units = "in", width = 5, height = 4, dpi = 300, compression = 'lzw')
}

# Assign the celltypes 

celltypes <- c("B", "CAF", "Macro", "NK", "T")

# Plot for each cell-type

for(cell_type in celltypes) {
  create_plots(cell_type, reference_x, X.matrix)
}

rodeo_data <- list("rodeo_matrix" = reference_x, "est_C_rodeo" = estimate_c_ref,
                   "cor_reference_x" = cor_cell_ref_test, "cor_X.matrix" = cor_cell_X_test)

saveRDS(rodeo_data, "Data/RODEO.rds")
