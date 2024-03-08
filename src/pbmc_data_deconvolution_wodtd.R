
# Load the pbmc data from *cor_estcwithDTD_pbmcpheno.R* file
pbmc_rodeo_ref_data <- readRDS(file.path("Data", "DTD_pbmc_rodeo_ref.rds"))

# Load the original pbmc data from the Data
pbmc_data <- readRDS("Data/pbmc_data.rds")

pbmc_train_expres <- pbmc_rodeo_ref_data$pbmc_train_expres
pbmc_test_expres <- pbmc_rodeo_ref_data$pbmc_test_expres
pbmc_train_pheno <- pbmc_data$train$bulk_pheno
pbmc_test_pheno <- pbmc_data$test$bulk_pheno
reference_profile_pbmc <- pbmc_rodeo_ref_data$reference_profile_pbmc

# Estimate cell-type proportion matrix using pbmc test data and estimated reference profile

library(DTD)


# Estimate cell proportion matrix without DTD

estimate_c_woDTD <- estimate_c(X.matrix = reference_profile_pbmc,
                                new.data = pbmc_test_expres,
                                DTD.model = rep(1, nrow(reference_profile_pbmc)),
                                estimate.c.type = "direct")

# Correlation between estimated c and pbmc test c

cor_estcwodtd_testcpbmc <- mapply(cor, as.data.frame(t(estimate_c_woDTD)),
                                  as.data.frame(t(pbmc_test_pheno)))
cor_estcwodtd_testcpbmc

# Make rds file 

rodeo_pbmc_data <- list( "naive_rodeo_reference" = cor_estcwodtd_testcpbmc)
                 


saveRDS(rodeo_pbmc_data, file.path("Data", "naive_pbmc_data.rds"))

