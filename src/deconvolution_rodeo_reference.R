

#### model deconvolution using rodeo generated reference matrix ####

# Load the rodeo data from *RODEO.R* file
rodeo_data <- readRDS("Data/RODEO.rds")

# Load the DTD data from *DTD tutorials.R* file 
dtd_data <- readRDS("Data/DTD tutorials.rds")

rodeo_reference <- rodeo_data$rodeo_matrix
training_data_dtd <- dtd_data$train
testing_data_dtd <- dtd_data$test


# Create a vector containing features from rodeo_reference

select_features <- rownames(rodeo_reference)

# Train deconvolution model using rodeo reference #

start_tweak <- rep(1, length(select_features))
names(start_tweak) <- select_features

# Deconvolute the model using X.matrix data 

library(DTD)

set.seed(100)

bulk_model <- train_deconvolution_model(
                                        tweak = start_tweak,
                                        X.matrix = rodeo_reference,
                                        train.data.list = training_data_dtd,
                                        test.data.list = training_data_dtd,
                                        estimate.c.type = "direct",
                                        use.implementation = "cpp")


# Estimate cell-type proportion matrix using new deconvolution model

estimate_c_rodeo_ref_model <- estimate_c(
                                        X.matrix = rodeo_reference,
                                        new.data = testing_data_dtd$mixtures,
                                        DTD.model = bulk_model,
                                        estimate.c.type = "direct")

# Correlation between estimated_c and DTD_train_pheno

cor_estc_dtd_testpheno <- mapply(cor, as.data.frame(t(testing_data_dtd$quantities)), 
                                  as.data.frame(t(estimate_c_rodeo_ref_model)))
cor_estc_dtd_testpheno

decon_rodeo_ref_data <- list("DTD_rodeo_reference" = cor_estc_dtd_testpheno)

saveRDS(decon_rodeo_ref_data, "Data/decon_rodeo_ref_data")
