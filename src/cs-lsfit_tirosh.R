
# Load the tirosh data from *DTD tutorial.R* file.

tiroshdata <- readRDS(file = "Data/DTD tutorials.rds")
training.data <- tiroshdata$train
testing.data <- tiroshdata$test
X.matrix <- tiroshdata$X_matrix

#### find reference matrix using expression and cellular composition matrix from naive training data ####

library(CellMix)

naive_lsfit_reference_data <- ged(training.data$mixtures, training.data$quantities, method = "cs-lsfit")

naive_lsfit_reference <- naive_lsfit_reference_data@fit@W

##### Estimating cell-composition matrix and find correlation without DTD using naive lsfit reference profile ####


#Load the DTD package

library(DTD)
set.seed(1234)

# Use test-expression matrix & estimated reference matrix from naive test data 

estimate_c_ref <- estimate_c(X.matrix = naive_lsfit_reference,
                             new.data = testing.data$mixtures,
                             DTD.model = rep(1,nrow(naive_lsfit_reference)),
                             estimate.c.type = "direct")

# Correlation between test and estimated_ref cell type composition matrices

cor_cell_ref_test <- mapply(cor, as.data.frame(t(estimate_c_ref)), 
                            as.data.frame(t(testing.data$quantities)))
cor_cell_ref_test



#### Estimating cell-composition matrix and find a correlation with DTD using naive lsfit reference profile ####

# Create a vector containing features from lsfit reference

select_features <- rownames(naive_lsfit_reference)

# Train deconvolution model using lsfit reference 

start_tweak <- rep(1, length(select_features))
names(start_tweak) <- select_features

# Deconvolute the model using naive lsfit reference matrix

library(DTD)

set.seed(100)

lsfit_bulk_model <- train_deconvolution_model(
                                              tweak = start_tweak,
                                              X.matrix = naive_lsfit_reference,
                                              train.data.list = training.data,
                                              test.data.list = training.data,
                                              estimate.c.type = "direct",
                                              use.implementation = "cpp")

# Estimate cell-type proportion matrix using new deconvolution model

estimate_c_lsfit_ref_model <- estimate_c(
                                          X.matrix = naive_lsfit_reference,
                                          new.data = testing.data$mixtures,
                                          DTD.model = lsfit_bulk_model,
                                          estimate.c.type = "direct")

# Correlation between estimated_c and DTD_test_pheno

cor_estc_dtd_testpheno <- mapply(cor, as.data.frame(t(testing.data$quantities)), 
                                 as.data.frame(t(estimate_c_lsfit_ref_model)))
cor_estc_dtd_testpheno



##### Find reference profile using bulk-expression data from Tirosh 19 tumors data #####

# Load filtered tirosh bulk expression data from *RODEO_Bulk_expres_phenodata.R* file

tirosh_bulk_data <- readRDS("Data/rodeo_bulkdata")
train_bulk_expres <- tirosh_bulk_data$rodeo_bulk_express
train_bulk_pheno <- tirosh_bulk_data$rodeo_bulk_pheno
test_bulk_expres <- tirosh_bulk_data$rodeo_test_bulk_express
test_bulk_pheno <- tirosh_bulk_data$rodeo_test_bulk_pheno

# find a reference matrix using cs-lsfit method

bulk_lsfit_reference_data <- ged(train_bulk_expres, t(train_bulk_pheno), method = "cs-lsfit")
bulk_lsfit_reference <- bulk_lsfit_reference_data@fit@W


#### Estimate cell-type proportion matrix without DTD using estimated bulk reference profile ####

estimatec_lsfit_testbulk <- estimate_c(X.matrix = bulk_lsfit_reference,
                                       new.data = test_bulk_expres,
                                       DTD.model = rep(1,nrow(bulk_lsfit_reference)),
                                       estimate.c.type = "direct")


#### Correlation between estimated C and bulk-pheno of test data ####

cor_estc_lsfit_testpheno <- mapply(cor, as.data.frame(t(estimatec_lsfit_testbulk)),
                                   as.data.frame(test_bulk_pheno))
cor_estc_lsfit_testpheno



#### Estimate cell-type proportion matrix with DTD using estimated bulk reference profile ####

# Generate bulk expression training data

bulk_training_data <- list("mixtures" = train_bulk_expres, "quantities" = t(train_bulk_pheno))

# Generate bulk model using bulk training data

bulk_select_features <- rownames(bulk_lsfit_reference)

bulk_start_tweak <- rep(1, length(bulk_select_features))
names(bulk_start_tweak) <- bulk_select_features

set.seed(100)

bulk_model_lsfitbulk <- train_deconvolution_model(
                                                    tweak = bulk_start_tweak,
                                                    X.matrix = bulk_lsfit_reference,
                                                    train.data.list = bulk_training_data,
                                                    test.data.list = bulk_training_data,
                                                    estimate.c.type = "direct",
                                                    use.implementation = "cpp"
                                                  )

# Estimate cell-type proportion matrix using new deconvolution model

estimate_c_bulk.model_lsfit <- estimate_c(
                                      X.matrix = bulk_lsfit_reference,
                                      new.data = test_bulk_expres,
                                      DTD.model = bulk_model_lsfitbulk$best.model$Tweak,
                                      estimate.c.type = "direct"
                                    )

# Find a correlation between estimated cell-type matrix and rodeo cell-type matrix

cor_estC_bulktest <- mapply(cor, as.data.frame(test_bulk_pheno),
                            as.data.frame(t(estimate_c_bulk.model_lsfit)))
cor_estC_bulktest


tirosh_lsfit_data <- list("naive_lsfit_reference" = cor_cell_ref_test, "DTD_lsfit_reference" = cor_estc_dtd_testpheno,
                          "naive_lsfit_reference_bulk" = cor_estc_lsfit_testpheno, "DTD_lsfit_reference_bulk" = cor_estC_bulktest)

saveRDS(tirosh_lsfit_data, "Data/tirosh_lsfit_data.rds")
