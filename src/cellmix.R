


####### downloading CellMix package #####


if (!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools", repos = c(CRAN = "https://cloud.r-project.org"))
}
if (!"csSAM" %in% installed.packages()[, "Package"]) {
  devtools::install_url("https://cran.r-project.org/src/contrib/Archive/csSAM/csSAM_1.2.4.tar.gz")
}
devtools::install_github("toscm/CellMix")


#####


X <- rmix(3, 100, 20, markers=5)
dim(X)
View(X)
