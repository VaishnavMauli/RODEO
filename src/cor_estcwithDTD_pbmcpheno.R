# Load the original pbmc data from Data

data_pbmc <- readRDS(file = "Data/pbmc_data.rds")

pbmc_train_expres <- data_pbmc$train$bulk_expr
pbmc_train_pheno <- data_pbmc$train$bulk_pheno
pbmc_test_pheno <- data_pbmc$test$bulk_pheno
pbmc_test_expres <- data_pbmc$test$bulk_expr



# E = expression matrix, C = cell type proportion matrix, X = reference matrix
library(MASS)
library(Rodeo)

# Estimate reference profile using pbmc train data 
reference_profile_pbmc <- Rodeo(E = pbmc_train_expres, C = pbmc_train_pheno)

# Create the number of features 

n.features <- 500

# Create vector with decreasing expression values
library(matrixStats)
sds_in_ref <- rowSds(reference_profile_pbmc)
names(sds_in_ref) <- rownames(reference_profile_pbmc)
sort_sds_in_ref <- sort(sds_in_ref,decreasing = TRUE)
select_features <- names(sort_sds_in_ref)[1:n.features]

##### Train a deconvolution model using rodeo bulk pbmc data ####

start_tweak <- rep(1, n.features)
names(start_tweak) <- select_features

library(DTD)

# Prepare pbmc training data using pbmc rodeo data

ref_rodeo_pbmc_top <- reference_profile_pbmc[select_features,]
pbmc_train_expres_top <- pbmc_train_expres[select_features,]
pbmc_test_expres_top <- pbmc_test_expres[select_features,]


training.data <- list("mixtures" = pbmc_train_expres_top, "quantities" = pbmc_train_pheno)

# Create model using pbmc training data as train as well as pbmc test data

set.seed(100)

pbmc_model <- train_deconvolution_model(
                                        tweak = start_tweak,
                                        X.matrix = ref_rodeo_pbmc_top,
                                        train.data.list = training.data,
                                        test.data.list = training.data,
                                        estimate.c.type = "direct",
                                        use.implementation = "cpp"
                                        )

# Estimate C using  pbmc data

est_c_withDTD <- estimate_c( X.matrix = ref_rodeo_pbmc_top,
                             new.data = pbmc_test_expres_top,
                             DTD.model = pbmc_model$best.model$Tweak,
                             estimate.c.type = "direct")

# correlation between estimated c and test c 

cor_estcwithDTD_testpbmc <- mapply(cor,
                                   as.data.frame(t(pbmc_test_pheno)),
                                   as.data.frame(t(est_c_withDTD)))
cor_estcwithDTD_testpbmc

# Create the data in a list and save as RDS file

pbmc_rodeo_data <- list("DTD_pbmc_rodeo_reference" = cor_estcwithDTD_testpbmc, "pbmc_train_expres" = pbmc_train_expres_top,
                        "reference_profile_pbmc" = ref_rodeo_pbmc_top, "pbmc_test_expres"= pbmc_test_expres_top)

saveRDS(pbmc_rodeo_data, file.path("Data", "DTD_pbmc_rodeo_ref.rds"))
