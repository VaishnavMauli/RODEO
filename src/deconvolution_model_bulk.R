

##### Create train deconvolution model using bulk expression training data ####

# Load the rodeo bulk data from *RODEO_Bulk_expres_phenodata.R* file

rodeo_bulkdata <- readRDS("Data/rodeo_bulkdata")
rodeo_bulk_expres <- rodeo_bulkdata$rodeo_bulk_express
rodeo_bulk_pheno <- rodeo_bulkdata$rodeo_bulk_pheno
rodeo_bulk_reference <- rodeo_bulkdata$reference_bulk
rodeo_test_bulk_expres <- rodeo_bulkdata$rodeo_test_bulk_express
rodeo_test_bulk_pheno <- rodeo_bulkdata$rodeo_test_bulk_pheno

# Create number of features 
n.features <- 500

# Create vector containing features with decreasing expression
library(matrixStats)

sds_in_ref <- rowSds(rodeo_bulk_reference)
names(sds_in_ref) <- rownames(rodeo_bulk_reference)
sort_sds_in_ref <- sort(sds_in_ref,decreasing = TRUE)
select_features <- names(sort_sds_in_ref)

#### Train a deconvolution model using bulk data ####

start_tweak <- rep(1, n.features)
names(start_tweak) <- select_features

library(DTD)

# Prepare a training data using rodeo data

training.data <- list("mixtures" = rodeo_bulk_expres, "quantities" = t(rodeo_bulk_pheno))

# Create model using training data as train as well as test data

set.seed(100)

bulk_model <- train_deconvolution_model(
              tweak = start_tweak,
              X.matrix = rodeo_bulk_reference,
              train.data.list = training.data,
              test.data.list = training.data,
              estimate.c.type = "direct",
              use.implementation = "cpp"
)

# Estimate cell-type proportion matrix using new deconvolution model

estimate_c_bulk.model <- estimate_c(
                          X.matrix = rodeo_bulk_reference,
                          new.data = rodeo_test_bulk_expres,
                          DTD.model = bulk_model,
                          estimate.c.type = "direct"
)

# Find a correlation between estimated cell-type matrix and rodeo cell-type matrix

cor_estC_bulktest <- mapply(cor, as.data.frame(rodeo_test_bulk_pheno),
                            as.data.frame(t(estimate_c_bulk.model)))
cor_estC_bulktest

decon_bulk_data <- list("cor_decon_model_bulk" = cor_estC_bulktest)

saveRDS(decon_bulk_data, "Data/decon_bulk_model")

