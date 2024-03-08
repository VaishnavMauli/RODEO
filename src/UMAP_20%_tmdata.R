
#### UMAP-part2 ####

# Load the data from DTD-tutorials from *DTD tutorials.R* file
tiroshdata <- readRDS(file = "Data/DTD tutorials.rds")

tm_pheno_data <- tiroshdata$tm_pheno_data
tm_expres <- tiroshdata$tm_expres

# Load the data from *RODEO.R* file
rodeodata <- readRDS(file = "Data/RODEO.rds")

rodeo_mat <- rodeodata$rodeo_matrix


set.seed(100)

#### 20% T cell-type samples ####

tm_pheno_t_pos <- which(tm_pheno_data$CellType  == "T")

tm_pheno_t_matrix <- as.matrix(tm_pheno_data[tm_pheno_t_pos,])


tm_pheno_t_data <- sample(x = 1:nrow(tm_pheno_t_matrix),
                          size = 0.2*nrow(tm_pheno_t_matrix))

t_mat <- as.matrix(tm_pheno_t_matrix[tm_pheno_t_data,])

#### 20% Unknown cell-type samples ####

tm_pheno_un_pos <- which(tm_pheno_data$CellType  == "unknown")

tm_pheno_un_matrix <- as.matrix(tm_pheno_data[tm_pheno_un_pos,])


tm_pheno_un_data <- sample(x = 1:nrow(tm_pheno_un_matrix),
                          size = 0.2*nrow(tm_pheno_un_matrix))

un_mat <- as.matrix(tm_pheno_un_matrix[tm_pheno_un_data,])

#### 20% NK cell-type samples ####

tm_pheno_nk_pos <- which(tm_pheno_data$CellType  == "NK")

tm_pheno_nk_matrix <- as.matrix(tm_pheno_data[tm_pheno_nk_pos,])

tm_pheno_nk_data <- sample(x = 1:nrow(tm_pheno_nk_matrix),
                          size = 0.2*nrow(tm_pheno_nk_matrix))

nk_mat <- as.matrix(tm_pheno_nk_matrix[tm_pheno_nk_data,])

#### 20% B cell-type samples ####

tm_pheno_b_pos <- which(tm_pheno_data$CellType  == "B")

tm_pheno_b_matrix <- as.matrix(tm_pheno_data[tm_pheno_b_pos,])

tm_pheno_b_data <- sample(x = 1:nrow(tm_pheno_b_matrix),
                          size = 0.2*nrow(tm_pheno_b_matrix))

b_mat <- as.matrix(tm_pheno_b_matrix[tm_pheno_b_data,])

#### 20% Macro cell-type samples ####

tm_pheno_m_pos <- which(tm_pheno_data$CellType  == "Macro")

tm_pheno_m_matrix <- as.matrix(tm_pheno_data[tm_pheno_m_pos,])

tm_pheno_m_data <- sample(x = 1:nrow(tm_pheno_m_matrix),
                          size = 0.2*nrow(tm_pheno_m_matrix))

m_mat <- as.matrix(tm_pheno_m_matrix[tm_pheno_m_data,])


# Merge all 20% cell-type samples in one matrix

new_pheno_mat <- rbind(t_mat, un_mat,nk_mat, b_mat, m_mat)


# Create new tm expression data with only 20% sample names

new_tm_com_data <- intersect(colnames(tm_expres), colnames(t(new_pheno_mat)))

new_tm_expres <- as.matrix(tm_expres[,new_tm_com_data])

new_tm_pheno <- as.matrix(new_pheno_mat[new_tm_com_data,])

# Find common genes from new tm expression data and Rodeo generated matrix

com_gene_data <- intersect(rownames(new_tm_expres), rownames(rodeo_mat))

filter_tm_expres <- as.matrix(new_tm_expres[com_gene_data,])


# Create new matrix by combining columns of filtered tm expres and rodeo matrix

matrix_tm_rodeo <- cbind(filter_tm_expres, rodeo_mat)


# Reduce dimentionality of combined matrix using UMAP function

library(uwot)
library(ggplot2)

set.seed(123)

umap_data <- umap(t(matrix_tm_rodeo),n_neighbors = 6, learning_rate = 0.6, 
                  init = "random", n_epochs = 1000, metric = "correlation")


# Create a data-frame to make a plot

tirosh_names <- rep("tirosh", 904)
rodeo_names <- rep("rodeo", 5)
group_data <- rbind(as.matrix(tirosh_names),as.matrix(rodeo_names))
celltype_data <- rbind(as.matrix(new_tm_pheno[,3]),as.matrix(colnames(rodeo_mat)))

data_plot <- data.frame(UMAP1  = umap_data[,1],
                        UMAP2 = umap_data[,2], 
                        Sample_name = colnames(matrix_tm_rodeo),
                        Celltype = celltype_data[,1],
                        Group = group_data[,1]
)
View(data_plot)                       

umap_plot <- ggplot(data = data_plot,
                    mapping = aes(x = UMAP1, y = UMAP2, group = Group))
umap_plot+
  geom_point(aes(color = Celltype, shape = Group), size =1.5) +
  scale_shape_manual(values = c(17,20)) +
  labs(x = "UMAP1", y ="UMAP2", subtitle = "UMAP_20%_tmdata")

 ggsave("UMAP_plot_tmdata.tiff", path = "Results", units="in", width=5, height=4, dpi=300, compression = 'lzw')


 tm_umap_data <- list("filter_tm_expres" = filter_tm_expres)
 saveRDS(tm_umap_data, "Data/tm_umap_data")
 
 