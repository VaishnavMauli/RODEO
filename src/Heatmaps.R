

##### Correlation heatmaps ####

#### Create heatmaps using $Tirosh data$ ####

#Load the correlation data

# Load data from *RODEO.R* file
rodeo_data <- readRDS("Data/RODEO.rds")

# Load data from *RODEO_Bulk_expres_phenodata.R* file
rodeo_bulk_data <- readRDS("Data/rodeo_bulkdata")

# Load data from *deconvolution_model_bulk.R* file
deconvoluion_bulk_data <- readRDS("Data/decon_bulk_model")

# Load data from *deconvlution_X.matrix.R* file
deconvolution_x.matrix_data <- readRDS("Data/decon_x.matrix_data")

# Load data from *deconvolution_rodeo_reference.R* file
decon_rodeo_ref_data <- readRDS("Data/decon_rodeo_ref_data")

# Load data from *cs-lsfit_tirosh.R* file
lsfit_data <- readRDS(file.path("Data", "tirosh_lsfit_data.rds"))

# Load data from *cs-qprog_tirosh.R* file
qprog_data <- readRDS(file.path("Data", "tirosh_qprog_data.rds"))


# Assign the correlation names

naive_rodeo_reference <- rodeo_data$cor_reference_x
naive_sc_reference <- rodeo_data$cor_X.matrix
naive_rodeo_reference_bulk <- rodeo_bulk_data$cor_reference_bulk
naive_lsfit_reference <- lsfit_data$naive_lsfit_reference
naive_qprog_reference <- qprog_data$naive_qprog_reference
naive_lsfit_reference_bulk <- lsfit_data$naive_lsfit_reference_bulk
naive_qprog_reference_bulk <- qprog_data$naive_qprog_reference_bulk
DTD_rodeo_reference_bulk <- deconvoluion_bulk_data$cor_decon_model_bulk
DTD_sc_reference <- deconvolution_x.matrix_data$DTD_sc_reference
DTD_rodeo_reference <- decon_rodeo_ref_data$DTD_rodeo_reference
DTD_lsfit_reference <- lsfit_data$DTD_lsfit_reference
DTD_qprog_reference <- qprog_data$DTD_qprog_reference
DTD_lsfit_reference_bulk <- lsfit_data$DTD_lsfit_reference_bulk
DTD_qprog_reference_bulk <- qprog_data$DTD_qprog_reference_bulk

# Prepare a matrix by combining columnwise to compare rodeo performance with lsfit and qprog algorithm
# Prepare the matrix without DTD and with DTD for rodeo, cs-lsfit and cs-qprog.

cor_data <- cbind(as.matrix(naive_rodeo_reference),
                  as.matrix(DTD_rodeo_reference),
                  as.matrix(naive_lsfit_reference),
                  as.matrix(DTD_lsfit_reference),
                  as.matrix(naive_qprog_reference),
                  as.matrix(DTD_qprog_reference),
                  as.matrix(naive_sc_reference),
                  as.matrix(DTD_sc_reference))
                    
# Give columnnames to cor_data matrix 
colnames(cor_data) <- c("naive_rodeo_reference", "DTD_rodeo_reference", "naive_lsfit_reference","DTD_lsfit_reference",
                        "naive_qprog_reference","DTD_qprog_reference","naive_sc_reference", "DTD_sc_reference")

# Prepare the matrix columnwise from bulk expression data of 19 tumors
cor_bulk_data <- cbind(as.matrix(naive_rodeo_reference_bulk),
                       as.matrix(DTD_rodeo_reference_bulk),
                       as.matrix(naive_lsfit_reference_bulk),
                       as.matrix(DTD_lsfit_reference_bulk),
                       as.matrix(naive_qprog_reference_bulk),
                       as.matrix(DTD_qprog_reference_bulk)
                       )

# Give the column names to the cor_bulk_data
colnames(cor_bulk_data) <- c("naive_rodeo_reference_bulk", "DTD_rodeo_reference_bulk",
                             "naive_lsfit_reference_bulk", "DTD_lsfit_reference_bulk",
                             "naive_qprog_reference_bulk", "DTD_qprog_reference_bulk")

# Load the required packages

library(reshape2)
library(ggplot2)

# Prepare a correlation matrix 

cor_matrix <- melt(cor_data)

# Create a heatmap using ggplot2

cor_plot <- ggplot(cor_matrix, mapping = aes(x = Var2, y = Var1, fill = value))

cor_plot +
  geom_tile(color = "white", size = 0.2) +
  geom_point(aes(size = abs(value)), shape = 15, show.legend = TRUE)+
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits = c(-1,1), name = "Correlation Value")+
  labs( subtitle = "Tirosh_Correlation_Heatmap", x = "Model names", y = " Cell types")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                   size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
  ggsave("heatmap1.tiff",
       path = "Results",
       units = "in", width = 6, height = 5, dpi = 300, compression = 'lzw')
  

# Deconvolution correlation matrix from bulk expression data

cor_bulk_matrix <- melt(cor_bulk_data)

# Create a heatmap using ggplot2

cor_bulk_plot <- ggplot(cor_bulk_matrix, mapping = aes(x = Var2, y = Var1, fill = value))

cor_bulk_plot +
  geom_tile(color = "white", size = 0.2) +
  geom_point(aes(size = abs(value)), shape = 15, show.legend = TRUE)+
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits = c(-1,1), name = "Correlation Value")+
  labs( subtitle = "Tirosh_Bulk_Correlation_heatmap", x = "Model names", y = " Cell types")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
  ggsave("heatmap2.tiff",
       path = "Results",
       units = "in", width = 6, height = 5, dpi = 300, compression = 'lzw')


  
#### Create heatmap using $PBMC data$ ####
  
# Load pbmc correlation data from *cs-qprog_pbmc.R* file
pbmc_qprog_data <- readRDS(file.path("Data", "pbmc_qprog_data.rds"))

# Load pbmc correlation data from *cs-lsfit_pbmc.R* file  
pbmc_lsfit_data <- readRDS(file.path("Data", "pbmc_lsfit_data.rds"))

# Load pbmc correlation data from *pbmc_data_deconvolution_wodtd.R* file
pbmc_rodeo_naive_data <- readRDS("Data/rodeo_pbmc_data")

# Load pbmc correlation data from *cor_estcwithDTD_pbmcpheno.R* file
pbmc_rodeo_dtd_data <- readRDS(file.path("Data","DTD_pbmc_rodeo_ref.rds"))

# Assign the data to variable

naive_pbmc_rodeo_reference <- pbmc_rodeo_naive_data$naive_rodeo_reference
naive_pbmc_lsfit_reference <- pbmc_lsfit_data$naive_lsfit_reference
naive_pbmc_qprog_reference <- pbmc_qprog_data$naive_qprog_reference
DTD_pbmc_rodeo_reference <- pbmc_rodeo_dtd_data$DTD_pbmc_rodeo_reference
DTD_pbmc_lsfit_reference <- pbmc_lsfit_data$DTD_lsfit_reference
DTD_pbmc_qprog_reference <- pbmc_qprog_data$DTD_qprog_reference

# Prepare a columnwise matrix with DTD and without DTD from rodeo, cs-lsfit and cs-qprog data

pbmc_cor_data <- cbind(as.matrix(naive_pbmc_rodeo_reference),
                       as.matrix(DTD_pbmc_rodeo_reference),
                       as.matrix(naive_pbmc_lsfit_reference),
                       as.matrix(DTD_pbmc_lsfit_reference),
                       as.matrix(naive_pbmc_qprog_reference),
                       as.matrix(DTD_pbmc_qprog_reference))

# Assign the name for columns of the generated matrix

colnames(pbmc_cor_data) <- c("naive_rodeo_reference", "DTD_pbmc_rodeo_reference",
                             "naive_lsfit_reference", "DTD_lsfit_reference",
                             "naive_qprog_reference", "DTD_qprog_reference")

# Deconvolution correlation matrix

pbmc_cor_matrix <- melt(pbmc_cor_data)

# Create a heatmap using ggplot2

pbmc_cor_plot <- ggplot(pbmc_cor_matrix, mapping = aes(x = Var2, y = Var1, fill = value))

pbmc_cor_plot +
  geom_tile(color = "white", size = 0.2) +
  geom_point(aes(size = abs(value)), shape = 15, show.legend = TRUE)+
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits = c(-1,1), name = "Correlation Value")+
  labs( subtitle = "PBMC_Correlation_Heatmap", x = "Model names", y = " Cell types")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
ggsave("heatmap3.tiff",
       path = "Results",
       units = "in", width = 6, height = 5, dpi = 300, compression = 'lzw')

