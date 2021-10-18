#' SRS
#' Self representation smoothing of single cell Hi-C matrix.
#'
#' @param scHiC The single-cell Hi-C matrix. It can take three types of formats. The preferred format is a single-cell matrix with
#' each column being a vector of the upper triangular matrix without including the diagonal entries
#' of the 2D matrix of a single-cell. Another types of formats are a list with each element being a
#' 2D single-cell contact matrix, or a 3D (\eqn{n\times n\times k}) array that has k matrices of dimension
#' \eqn{n\times n}. scHiCSRS automatically transforms these two types of input into a matrix with each
#' column being the vector of upper triangular matrix of a single-cell. For a single-cell matrix of
#' size \eqn{n \times n}, the length of the vector should be \eqn{n\times(n-1)/2}. We only need the upper
#' triangular matrix because the Hi-C matrix are symmetrical.
#' @param expected Underline true counts of the simulated data. For real data analysis, just set it as NULL.It takes
#'three formats that is the same as scHiC.
#' @param windowsize The size of neighborhood region. A windowsize of w results in a (2w+1)*(2w+1) neighboring submatirx.
#' @param nbins Number of bins of the observed single cell HiC matrix.
#' @param lambda1 Tuning parameter to facilitate feature selection and regularization.
#' @param lambda2 Tuning parameter to penalize teh diagonal element of the parameter to eliminate the trivial solution of representing an expression level as a linear combination of itself.
#' @param initA The initialization of A. The elements of A represent the similarities between loci in the same cell.
#' @param initS The initialization of S. The elements of S represent the similarities between all single cells at the same position.
#' @param ncores Number of cores to use. Default is 1.
#' @param MAX_ITER Maximum iteration of the external circulation of SRS.
#' @param ABSTOL Absolute tolerance of the external circulation.
#' @param learning_rate A hyper parameter that controls the speed of adjusting the weights of the network with respect to the loss gradient.
#' @param epochs The number of the entire training set going through the entire network.
#' @param batch_size The number of examples that are fed to the algorithm at a time.
#' @param run_batch Whether to use batch or to set the number of all the samples as teh value of the batch size. Default is TRUE.
#' @param verbose Whether to output the value of metrics at the end of each epoch. Default is TRUE.
#' @param estimates.only If TRUE, than out the SRS imputed matrix. If FALSE, A list of information is outputted.
#'
#'
#' @return
#' @export
#'
#' @import SAVER
#'
#' @import keras
#'
#' @import tensorflow
#'
#' @import mclust
#'
#' @import reticulate
#'
#' @examples
#' SRS(scHiC, windowsize=2, nbins=61, learning_rate = 0.0001,epochs = 100)

SRS <- function(scHiC, expected, windowsize=2, nbins, lambda1 = NULL, lambda2 = 1e10, initA = NULL,

                initS = NULL, ncores = 1, MAX_ITER = 4, ABSTOL = 1e-3, learning_rate = 0.0001,

                  epochs = 100, batch_size = 128, run_batch = TRUE, verbose = TRUE, estimates.only = FALSE){

  ##readin data of different types
  if(is.list(scHiC)){
    single=NULL
    for (c in 1:length(scHiC)) {
      single=cbind(single, mattovec(scHiC[[c]]))
    }
  }else if(is.array(scHiC) & length(dim(scHiC))>2){
    single=NULL
    for (c in 1:(dim(scHiC)[3])) {
      single=cbind(single, mattovec(scHiC[,,c]))
    }
  }else if(is.matrix(scHiC)){
    single=scHiC
  }else{
    stop("The input type must be a matrix, 3D array or a list of 2D matrix.")
  }

  if (nrow(single) != (nbins*(nbins-1)/2)){
    stop("The provided n and single cells do not match in size!")
  }

  scHiC=single

  # calculating the index matrix
  nei_index=NULL
  for (j in 2:nbins) {
    for (i in 1:(j-1)) {
      mm=matrix(0,nbins,nbins)
      mm[max(1,i-windowsize):min(nbins,i+windowsize),max(1,j-windowsize):min(nbins,j+windowsize)]=1
      nei_index=cbind(nei_index,mm[upper.tri(mm, diag = FALSE)])
    }
  }
  nei_index=(nei_index==0)
  nei_index = t(nei_index)
  # calculating the tunning parameter lambda

  if (is.null(lambda1)) {

    #message("Calculating the penalty parameter lambda ...")

    lambda1 <- calc.lambda(scHiC)

    #message("Done!")

  }


  # preprocessing and log-normalization

  #message("Starting preprocessing and log-normalization ...")

  X <- log_normalization(scHiC)

  #message("Done!")


  scHiC <- clean.data(scHiC)

  npairs <- nrow(scHiC)

  ncells <- ncol(scHiC)

  gene.names <- rownames(scHiC)

  cell.names <- colnames(scHiC)


  # assign size factor

  sf.out <- calc.size.factor(scHiC)

  sf <- sf.out[[1]]

  scale.sf <- sf.out[[2]]



  result <- list()

  result$Impute_All <- matrix(0, npairs, ncells, dimnames = list(gene.names, cell.names))

  result$A <- matrix(0, npairs, npairs)

  result$S <- matrix(0, ncells, ncells)

  #result$obj.prior <- c()  #############



  # if (!estimates.only) {
  #
  #   result$se <- matrix(0, npairs, ncells, dimnames = list(gene.names, cell.names))
  #
  # } else {
  #
  #   result$se <- NA
  #
  # }

  #result$info <- c(list(0), list(0), list(0), list(0))

  #names(result$info) <- c("size.factor", "pred.time", "posterior.time", "total.time")

  #result$info$size.factor <- scale.sf*sf



  # initialize A and S

  if (is.null(initA)) {

    initA <- matrix(0, npairs, npairs)

  }

  if (is.null(initS)) {

    initS <- matrix(0, ncells, ncells)

  }

  A <- initA

  S <- initS


  #message("Imputation starts ...")

  #message("Calculating the prior mean for each gene in each cell ...")

  #message("iter", ' objective')

  pred.st <- Sys.time()

  k <- 1

  history.objval <- c()

  zeros = matrix(0, npairs, ncells)

  while (k <= MAX_ITER){

    Z <- pmax(X, zeros)

    Y <- X - A%*%X

    S_old <- S

    if (!run_batch){

      batch_size = nrow(Z)

    } else{

      batch_size = batch_size

    }

    # using keras with SGD algorithm to calculate S or A

    S <- keras_lasso_regression(Y, Z, epochs = epochs, batch_size = batch_size, lambda1 = lambda1, lambda2 = lambda2,

                                learning_rate = learning_rate, verbose = verbose)

    if (k > 1){

      S <- (S_old + S)/2

    } else {

      S <- S

    }

    Z <- t(pmax(X, zeros))

    Y <- X - X%*%S

    A_old <- A

    if (!run_batch){

      batch_size = nrow(Z)

    } else{

      batch_size = batch_size

    }


    A <- t(keras_lasso_regression(t(Y), Z, epochs = epochs, batch_size = batch_size, lambda1 = lambda1, lambda2 = lambda2,

                                  learning_rate = learning_rate, verbose = verbose))


    if (k > 1){

      A <- (A_old + A)/2

    } else {

      A <- A

    }

    A[nei_index]=0 #force non-neighbors to be 0

    history.objval[k] <- objective(X, lambda1, A, S)

    #message(k, '    ', round(history.objval[k], 3))

    if (k > 1 && abs(history.objval[k] - history.objval[k-1]) < ABSTOL) break

    k <- k + 1

  }

  Xhat <- pmax(A%*%X + X%*%S, zeros)

  obj.prior <- history.objval #########


  pred.time <- Sys.time() - pred.st

  message("Calculating the posterior means with ", ncores, " worker(s)")

  posterior.st <- Sys.time()

  out <- SAVER::saver(scHiC, ncores = ncores, size.factor=sf, mu = exp(as.matrix(Xhat)))##investigate

  posterior.time <- Sys.time() - posterior.st

  total.time <- Sys.time() - pred.st

  estimate<- out$estimate%*%diag(sf)

  result$scHiC <- as.matrix(scHiC)

  result$Impute_All <- estimate

  #result$se <- out$se

  #result$info$pred.time <- pred.time

  #result$info$posterior.time <- posterior.time

  result$total.time <- total.time

  result$A <- A

  result$S <- S

  result$Expected <- expected

  SZ <- orga_mix(observed = scHiC, imputed = estimate)

  result$SZ <- SZ

  Impute_SZ <- estimate
  Impute_SZ[SZ>0.5] <- 0
  result$Impute_SZ <- Impute_SZ

  #result$obj_prior <- obj.prior    ###########

  #message("Done!")

  #message("Finish time: ", Sys.time())

  #message("Total time: ", format(result$info$total.time))

  if (!estimates.only) {

    class(result) <- "scHiCSRS"

    result

  } else {

    result$Impute_All

  }

}
