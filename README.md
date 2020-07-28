# **EV-origin**

#### **A computational approach to robustly enumerating the abundance of tissue-cellular origins using circulating exLR-seq profiles**


##### **Summary**

This repository contains the source R Codes and two reference signature matrixes of EV-origin program to respectively estimate the relative and absolute abundance of hemopoietic and tissue components from exLR-seq data, it also applies the detail description of model constructing process so that others can replicate our deconvolution results and use our candidate machine learning method on other types of EVs RNA sequencing data.

##### **Download and Install**

Software requirements: R version 3.4.2 or later (dependencies below might not work properly with earlier versions). It is recommended to manipulate our R script on windows platform. You can type the following command to download the latest R packages for EV-origin as follow:

install.packages("MASS")

install.packages("glmnet")

install.packages("pracma")

install.packages("nnls")

install.packages("e1071")

install.packages("parallel")

install.packages("preprocessCore")

After installing, you can respectively load these packages with R command library("xxx").

For more detailed usage instructions, see our manuscript on the journal of Molecules.

##### **Main strategy and detail description of model construction**

The concept of EV-origin deconvolution is to find the optimal solution of a convoluting equation expressed as AX = B, where A is the transcriptome mixture of the exLR-seq profile, B is the comparable signature matrix for the expression of genes in all types of subset, and X is the vector of relative/absolute proportions of all cell/tissue components. 

Additionally, with a final filtered nu-support vector regression (ν-SVR) model, our goal was to investigate a hyperplane that fits as many data points as possible within an optimal distance. Three main steps were included in our EV-origin process. The first step was the zero-mean normalization of the input exLR-seq expression data. The second step was parameter selection. The ν -SVR model with a linear kernel was tested with different values of ν (ranging from 0 to 1). The parameter with the lowest root mean square error (RMSE) was kept for variable shrinking and model construction. Finally, the relative/absolute proportions of tissue-cellular components in each sample were calculated based on these optimized parameters.

 
