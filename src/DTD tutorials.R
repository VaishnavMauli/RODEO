

### DTD Tutorial ###

#This tutorial is written by Marian Schoen on deconvolution model.



#Set the deconvolution r-package in local repository
library(DTD)

# for downloading model sc-RNA seq data from NCBI
library(GEOquery)

# set package to check on run-time of optimization
library(tictoc)

#Obtained the gene profile data from GEO (Gene Expression Omnibus.
#GEO is maintained by NCBI. 
rawdata <- getGEOSuppFiles( GEO = "GSE72056", makeDirectory = TRUE)

#The getGEOSuppFiles function creates a directory named "GSE72056", where it stored in .txt.gz. Read it in via:




tirosh.melanoma <- read.table(file = "Data/GSE72056/GSE72056_melanoma_single_cell_revised_v2.txt.gz", 
                               stringsAsFactors = F,
                               header = T,
                               sep = "\t")

#First 3 row holds
tm.pheno <- as.matrix(tirosh.melanoma[1:3, -1])
rownames(tm.pheno) <- tirosh.melanoma[1:3, 1]

# for duplicated rownames

# to convert the gene names into character
row.names <- as.character(tirosh.melanoma[4:nrow(tirosh.melanoma), 1])

# To find out at which positions duplicated genes are in row.names
dupli.pos <- which(duplicated(row.names))

# To concatenated duplicates genes into unique.names
unique.names <- paste0(row.names[dupli.pos], "--2")

row.names[dupli.pos] <- unique.names
 
# Undo log transformation from tirosh.melanoma (DTT is additive) 
# Our data is log2 based, so used 2^(data)-1

tm.expr <- as.matrix(2^(tirosh.melanoma[4:nrow(tirosh.melanoma), -1]) - 1)

#Reset the rownames

rownames(tm.expr) <- row.names

#Remove NA values from expression data

tm.expr[!is.finite(tm.expr)] <- 0



#Normalize each profile to fixed number of counts

tm.expr <- normalize_to_count(tm.expr)




## Phenotype information ##

# Given raw data is present has malignant and cell type rows 
# Which present in numeric data. We need to convert into strings.
#To convert malignant data into strings.

map.malignant <- function(x) {
  if (x ==1) return("Non-malignant")
  if (x ==2) return("malignant")
  if (x ==3) return("unspecified")
  return("unresolved")
}

#To convert cell-type data into strings.

map.celltype <- function(x) {
  if (x ==1) return("T")
  if (x ==2) return("B")
  if (x ==3) return("Macro")
  if (x ==4) return("Endo")
  if (x ==5) return("NK")
  if (x ==6) return("CAF")
  
  return("unknown")
}

# Make a data frame for a string

tm.pheno.readable <- data.frame(
  "tumor" = tm.pheno["tumor",],
  "malignanttype" = sapply(
    tm.pheno["malignant(1=no,2=yes,0=unresolved)",],
    map.malignant
  ),
  "CellType" = sapply(
    tm.pheno["non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)",],
    map.celltype
                     )
)
head(tm.pheno.readable)

# This is our data-matrix and now we don't need tm.pheno

rm(tm.pheno)

### Single cell profiles ####

head(apply(tm.expr, 2, sum))

# Remember in tm.expr, each row is feature and column is single cell profile

tm.expr[1:5, 1:2]


##### Reconstruct the inferred bulk profiles ####


tumor.names <- as.character(unique(tm.pheno.readable$tumor))

length(tumor.names)

# Initiate empty expression matrix #

bulk.expres <- matrix(NA, 
                      nrow = nrow(tm.expr), 
                      ncol = length(tumor.names))

# Naming the rows and columns of bulk.expres

rownames(bulk.expres) <- rownames(tm.expr)
colnames(bulk.expres) <- tumor.names

#Phenodata is obtained by given following code.
 
bulk.pheno <- matrix(0, 
                     nrow = length(tumor.names),
                     ncol = length(unique(tm.pheno.readable$CellType)))
rownames(bulk.pheno) <- tumor.names
colnames(bulk.pheno) <- unique(tm.pheno.readable$CellType)

# Iterate over each tumor and sum up all profiles

for ( L.tumor in tumor.names ){
  tmp.samples <- which(tm.pheno.readable$tumor == L.tumor)
  bulk.expres[, L.tumor] <- rowSums(tm.expr[, tmp.samples])
  
  tmp.table <- table(tm.pheno.readable[tmp.samples, "CellType"])
  bulk.pheno[L.tumor, names(tmp.table)] <- tmp.table/ sum(tmp.table)
}

# Normalize the profiles 

library(DTD)
bulk.expres <- normalize_to_count(bulk.expres)

##### Split scRNA seq data into training and test dataset ####

set.seed(100)
train.data <- sample(x = 1:length(tumor.names),
                     size = 0.5* length(tumor.names))

test.data <- sample(1:length(tumor.names)) [-train.data]

# Concatanate the name 

cat("Training tumors:", tumor.names[train.data])

cat("test tumors:", tumor.names[test.data] )

## Testing ##
train.melanomas <- tumor.names[train.data]
test.melanomas <- tumor.names[test.data]

#train data
train.profile.pos <- which(tm.pheno.readable$tumor %in% train.melanomas)
train.pheno <- tm.pheno.readable[train.profile.pos,]  
train.profiles <- tm.expr[, train.profile.pos]


#test data
test.profiles.pos <- which(tm.pheno.readable$tumor %in% test.melanomas)
test.pheno <- tm.pheno.readable[test.profiles.pos, ]
test.profiles <- tm.expr[, test.profiles.pos]



#### DTD Analysis ####

# Starting DTD analysis with creating vector
# Vector that maps profile to celltype

indicator.vector <- as.character(tm.pheno.readable$CellType)
names(indicator.vector) <- rownames(tm.pheno.readable)
indicator.train <- indicator.vector[train.profile.pos]

indicator.test <- indicator.vector[test.profiles.pos]


include.in.X <- c("B", "CAF", "Macro", "NK", "T")


## Generate a reference matrix X ##

sample.X <- sample_random_X(included.in.X = include.in.X,
                            pheno = indicator.train,
                            expr.data = train.profiles,
                            percentage.of.all.cells = 0.2, #percentage to generate cell type profile
                            normalize.to.count = TRUE
                           )

X.matrix <- sample.X$X.matrix  #function returns X.matrix


#Already used profiles must not use further
samples.to.remove <- sample.X$samples.to.remove


#They must be remove train.profiles
remain.train.profiles <- train.profiles[,
                                         -which(colnames(train.profiles) %in% samples.to.remove)]


#And corresponding indicator vector will be 
remain.indicator.train <- indicator.train[colnames(remain.train.profiles)]


# Now, we reduce the number of feature
n.features <- 500

#Pre-selection is done using s.d. of refernce matrix X
sds.in.x <- rowSds(X.matrix)
names(sds.in.x) <- rownames(X.matrix)
sorted.sds.in.x <- sort(sds.in.x, decreasing = TRUE)


#Select the top n.features 
select.features <- names(sorted.sds.in.x)[1:n.features] 

#Now, reduce the expression set for all matrices 
X.matrix <- X.matrix[select.features,]
remaining.train.profiles <- remain.train.profiles[select.features,]
test.profiles <- test.profiles[select.features,]
bulk.exprs <- bulk.expres[select.features,]


##### Generate training and test "in-silico" mixtures ####

n.samples <- n.features     

# There are nearly 2500 SC profiles in the training set of remain profile
# Choose n.per.mixture around ~20% of training set 
n.per.mixture <- 400

#Generate train mixtures

training.data <- mix_samples(expr.data = remaining.train.profiles,
                             pheno = remain.indicator.train,
                             included.in.X = include.in.X,
                             n.samples = n.samples,
                             n.per.mixture = 400,
                             verbose = FALSE
                             )
#generate test mixtures 
testing.data <- mix_samples(expr.data = test.profiles,
                             pheno = indicator.test,
                             included.in.X = include.in.X,
                             n.samples = n.samples,
                             n.per.mixture = n.per.mixture,
                             verbose = FALSE
)

tirosh_data <- list("test" = testing.data, "train" = training.data,
                    "bulk_expres" = bulk.exprs, "X_matrix" = X.matrix, "bulk_pheno" = bulk.pheno, 
                    "tumor_names"= tumor.names, "tm_pheno_data" = tm.pheno.readable, "tm_expres" = tm.expr)

saveRDS(tirosh_data,"Data/DTD tutorials.rds")
