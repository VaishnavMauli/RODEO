
# Load the pbmc data and assign the matrices

pbmc_data <- readRDS(file.path("Data", "pbmc_data.rds"))
pbmc_train <- pbmc_data$train
pbmc_test <- pbmc_data$test
pbmc_train_expres <- pbmc_train$bulk_expr
pbmc_train_pheno <- pbmc_train$bulk_pheno
pbmc_test_expres <- pbmc_test$bulk_expr
pbmc_test_pheno <- pbmc_test$bulk_pheno

# Load the CellMix package
library(CellMix)

# Estimate gene expression reference profile using ged function by cs-lsfit method

reference_lsfit_data <- ged(pbmc_train_expres, pbmc_train_pheno, method = "cs-lsfit")

reference_lsfit <- reference_lsfit_data@fit@W




#### Estimate cell-type proportion matrix with DTD and find a correlation ####

# Create the number of features 

n.features <- 500

# Create vector with decreasing expression values

library(matrixStats)
sds_in_ref <- rowSds(reference_lsfit)
names(sds_in_ref) <- rownames(reference_lsfit)
sort_sds_in_ref <- sort(sds_in_ref,decreasing = TRUE)
select_features <- names(sort_sds_in_ref)[1:n.features]

##### Train a deconvolution model using rodeo bulk pbmc data ####

start_tweak <- rep(1, n.features)
names(start_tweak) <- select_features

library(DTD)

# Prepare pbmc training data using pbmc rodeo data

ref_lsfit_pbmc_top <- reference_lsfit[select_features,]
pbmc_train_expres_top <- pbmc_train_expres[select_features,]
pbmc_test_expres_top <- pbmc_test_expres[select_features,]


training.data <- list("mixtures" = pbmc_train_expres_top, "quantities" = pbmc_train_pheno)

# Create model using pbmc training data 

set.seed(100)

pbmc_model <- train_deconvolution_model(
                                          tweak = start_tweak,
                                          X.matrix = ref_lsfit_pbmc_top,
                                          train.data.list = training.data,
                                          test.data.list = training.data,
                                          estimate.c.type = "direct",
                                          use.implementation = "cpp"
                                        )

# Estimate C using  pbmc data

est_c_withDTD <- estimate_c( X.matrix = ref_lsfit_pbmc_top,
                             new.data = pbmc_test_expres_top,
                             DTD.model = pbmc_model$best.model$Tweak,
                             estimate.c.type = "direct")

# correlation between estimated c and test c 

cor_estcwithDTD_testpbmc <- mapply(cor,
                                   as.data.frame(t(pbmc_test_pheno)),
                                   as.data.frame(t(est_c_withDTD)))
cor_estcwithDTD_testpbmc


###### Estimate cell-type proportion matrix without DTD and find a correlation #### 

# Use pbmc test data and estimated reference profile

library(DTD)


# Estimate cell proportion matrix without DTD

estimate_c_woDTD <- estimate_c(X.matrix = ref_lsfit_pbmc_top,
                               new.data = pbmc_test_expres_top,
                               DTD.model = rep(1, nrow(ref_lsfit_pbmc_top)),
                               estimate.c.type = "direct")

# Correlation between estimated c and pbmc test c

cor_estcwodtd_testcpbmc <- mapply(cor, as.data.frame(t(estimate_c_woDTD)),
                                  as.data.frame(t(pbmc_test_pheno)))
cor_estcwodtd_testcpbmc


pbmc_lsfit_data <- list("naive_lsfit_reference" = cor_estcwodtd_testcpbmc, "DTD_lsfit_reference" = cor_estcwithDTD_testpbmc)

saveRDS(pbmc_lsfit_data, "Data/pbmc_lsfit_data.rds")
