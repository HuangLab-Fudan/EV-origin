library(e1071)
### nu-SVR
### read the reference matrix as "ref.m"
x <- data.matrix(exp)
x <- (x - mean(x)) / sd(as.vector(x))
avdata.m <- x#the expression profile normalization
nu.v = c(0.25, 0.5, 0.75)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do.exLR-origin <- function(avdata.m, ref.m, nu.v) {
  map.idx <- match(rownames(ref.m), rownames(avdata.m))
  rep.idx <- which(is.na(map.idx) == FALSE)
  data2.m <- avdata.m[map.idx[rep.idx], ]
  ref2.m <- ref.m[rep.idx, ]
  est.ab.lm <- list()
  est.lm <- list()
  nui <- 1
  for (nu in nu.v) {
    est.m <- matrix(NA, nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.m) <- colnames(ref2.m)
    rownames(est.m) <- colnames(data2.m)
    est.ab.m <- matrix(NA, nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.ab.m) <- colnames(ref2.m)
    rownames(est.ab.m) <- colnames(data2.m)
    
    for (s in seq_len(ncol(data2.m))) {
      svm.o <- svm(x = ref2.m, y = data2.m[, s], scale = TRUE, type = "nu-regression", kernel = "linear", nu = nu)
      coef.v <- t(svm.o$coefs) %*% svm.o$SV
      coef.v[which(coef.v < 0)] <- 1*10^-10
      est.ab.m[s,] <- coef.v
      total <- sum(coef.v)
      coef.v <- coef.v/total
      est.m[s, ] <- coef.v
    }
    est.lm[[nui]] <- est.m
    est.ab.lm[[nui]] <- est.ab.m
    nui <- nui + 1
  }
  
  #### select best nu using RMSE
  rmse.m <- matrix(NA, nrow = ncol(avdata.m), ncol = length(nu.v))
  for (nui in seq_along(nu.v)) {
    reconst.m <- ref2.m %*% t(est.lm[[nui]])
    s <- seq_len(ncol(avdata.m))
    rmse.m[s, nui] <- sqrt(colMeans((data2.m[, s] - reconst.m[, s])^2))
    message(nui)
  }
  colnames(rmse.m) <- nu.v
  nu.idx <- apply(rmse.m, 1, which.min)
  estF.m <- est.m
  for (s in seq_len(nrow(estF.m))) {
    estF.m[s, ] <- est.lm[[nu.idx[s]]][s, ]
  }
  estF.ab.m <- est.ab.m
  for (s in seq_len(nrow(estF.ab.m))) {
    estF.ab.m[s, ] <- est.ab.lm[[nu.idx[s]]][s, ]
  }
  #selecting min RMSE
  rmse.min.value <- as.data.frame(apply(rmse.m, 1, min))
  rownames(rmse.min.value) <- colnames(avdata.m)
  colnames(rmse.min.value)[1] <- 'RMSE'
  #caculating PCC
  pearson.corr.value <- c()
  for (i in 1:ncol(data2.m)) {
    cor.index <- cor.test(data2.m[, i], reconst.m[, i])
    cor.p <- as.numeric(cor.index$estimate)
    pearson.corr.value <- c(pearson.corr.value, cor.p)
  }
  estF.m <- cbind.data.frame(estF.m, pearson.corr.value, rmse.min.value)
  return(list(estF = estF.m, est.ab.sum = estF.ab.m, nu = nu.v[nu.idx]))
}
exLR-origin.results <- do.exLR-origin(avdata.m, ref.m, nu.v)
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
