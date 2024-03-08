
# Load the data from *DTD tutorials.R* file 
tiroshdata <- readRDS(file = "Data/DTD tutorials.rds")

# Bulk expression and phenotype data 
tm_pheno_data <- tiroshdata$tm_pheno_data
Bulk_expres <- tiroshdata$bulk_expres
bulk.phenodata <- tiroshdata$bulk_pheno

#### Divide the bulk expression data in training and test-dataset ####

set.seed(12879)

tumor.names <- colnames(Bulk_expres)

train_bulk <- sample(x = 1:length(tumor.names),
                            size = 0.6*length(tumor.names))

test_bulk <- sample(x = length(tumor.names)) [-train_bulk]

#Put expression values in respective train and test dataset in a matrix                  
matrix_train_bulk_expres <- as.matrix(Bulk_expres[,train_bulk])

matrix_test_bulk_expres <- as.matrix(Bulk_expres[,test_bulk])


##### Divide the bulk phenotype data in training and test-dataset #####

matrix_train_bulk_pheno <- as.matrix(bulk.phenodata[train_bulk,])

matrix_test_bulk_pheno <- as.matrix(bulk.phenodata[test_bulk,])



#### Estimate the reference profile using training bulk expression data ####

library(MASS)
library(Rodeo)

reference_profile_bulk <- Rodeo(E = matrix_train_bulk_expres, C = t(matrix_train_bulk_pheno))
View(reference_profile_bulk)


#### Estimate Cell-type proportion matrix using test-bulk expression data ####

library(DTD)

estimateC_testbulk <- estimate_c(X.matrix = reference_profile_bulk,
                                 new.data = matrix_test_bulk_expres,
                                 DTD.model = rep(1,nrow(reference_profile_bulk)),
                                 estimate.c.type = "direct")


#### Correlation between estimated C and bulk-pheno of test data ####

cor_estimateC_testpheno <- mapply(cor, as.data.frame(t(estimateC_testbulk)),
                                       as.data.frame(matrix_test_bulk_pheno))
cor_estimateC_testpheno


### create RDS file to save the data ###

rodeo_bulkdata <- list("reference_bulk" =  reference_profile_bulk, 
                       "rodeo_bulk_express" = matrix_train_bulk_expres,
                       "rodeo_bulk_pheno" = matrix_train_bulk_pheno,
                       "rodeo_test_bulk_express" = matrix_test_bulk_expres,
                       "rodeo_test_bulk_pheno" = matrix_test_bulk_pheno,
                       "cor_reference_bulk" = cor_estimateC_testpheno)
saveRDS(rodeo_bulkdata, "Data/rodeo_bulkdata")
