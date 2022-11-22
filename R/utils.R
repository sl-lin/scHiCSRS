#######################################################
##              Gaussian mixture model               ##
#######################################################
### estimate parameters in the Gaussian mixture distribution

get_mix <- function(xdata) {

  bic=NULL
  for (g in 1:10) {
    fit = Mclust(xdata, G=g, model="V")
    bic[g] = fit$BIC
  }
  nmix=which.max(bic)

  fit = Mclust(xdata, G=nmix, model="V")
  paramt = c(fit$parameters$pro, fit$parameters$mean,fit$parameters$variance$sigmasq)
  names(paramt)=NULL
  return(paramt)
}


dmix <- function(x,pars) {
  output=c()
  G=length(pars)/3
  if(G==1){ ### G=1
    output=rep(1, length(x))
    return(output)
  }else if(G==2){ ### G=2
    for (j in 1:length(x)) {
      comp=c()
      for (i in 1:G) {
        comp[i]=pars[i]*dnorm(x[j], mean = pars[i+G], sd=sqrt(pars[i+2*G]))
      }
      output[j]=sum(comp[1])/sum(comp)
    }
    return(output)
  }else{ ### G>=3
    d=rep(0,(G-1))
    for (gg in 1:(G-1)) {
      d[gg]=pars[G+gg+1]-pars[G+gg]
    }

    g=1
    ind=NULL
    while (g < G) {
      if(10*d[g] > d[g+1]){
        ind=c(ind,g)
        break
      }else{
        ind=c(ind, g, g+1)
      }
      g=g+1
    }

    for (j in 1:length(x)) {
      comp=c()
      for (i in 1:G) {
        comp[i]=pars[i]*dnorm(x[j], mean = pars[i+G], sd=sqrt(pars[i+2*G]))
      }
      output[j]=sum(comp[ind])/sum(comp)
    }
    return(output)
  }#end of else
}


orga_mix <- function(observed, imputed) {

    imputed = lnorm(imputed)

    prob = matrix(0, nrow = nrow(imputed), ncol = ncol(imputed))
    prob[which(apply(observed,1,sum)==0),]=1
    dat1=observed
    dat1[which(apply(observed,1,sum)==0),]=1
    index=(dat1==0)
    dat2=imputed[index]
    paramt=get_mix(dat2)
    pp=dmix(dat2,paramt)
    prob[index]=pp

    return(prob)
}



# dmix <- function(x,pars) {
#   output=c()
#   G=length(pars)/3
#   d=rep(0,(G-1))
#   for (gg in 1:(G-1)) {
#     d[gg]=pars[G+gg+1]-pars[G+gg]
#   }
#   if((10*d[1])<=d[2]){
#     ind=c(1,2)
#   }else{
#     if((10*d[2])<=d[3]){
#       ind=c(1,2,3)
#     }else{
#       ind=c(1)
#     }
#    }
#   for (j in 1:length(x)) {
#     comp=c()
#     for (i in 1:G) {
#       comp[i]=pars[i]*dnorm(x[j], mean = pars[i+G], sd=sqrt(pars[i+2*G]))
#     }
#     output[j]=sum(comp[ind])/sum(comp)
#   }
#   return(output)
# }


get_mix_parameters <-
  function (count, nmix, ncores = 1)
  {
    count = as.matrix(count)
    #null_genes = which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
    parslist = mclapply(1:nrow(count), function(ii) {
      #  if (ii %% 2000 == 0) {
      #    gc()
      #    print(ii)
      #  }
      #  if (ii %in% null_genes) {
      #   return(rep(NA, 5))
      #  }
      xdata = count[ii, ]
      paramt = try(get_mix(xdata, nmix), silent = TRUE)
      if (class(paramt) == "try-error"){
        paramt = rep(NA, 5)
      }
      return(paramt)
    }, mc.cores = ncores)
    parslist = Reduce(rbind, parslist)
    #colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
    return(parslist)
  }


lnorm <- function(raw_count) {
  totalCounts_by_cell = colSums(raw_count)
  totalCounts_by_cell[totalCounts_by_cell == 0] = 1
  raw_count = sweep(raw_count, MARGIN = 2, median(totalCounts_by_cell)/totalCounts_by_cell, FUN = "*")
  if (min(raw_count) < 0) {
    stop("smallest read count cannot be negative!")
  }
  count_lnorm = log10(raw_count + 1.01)
  return(count_lnorm)
}

matrow <- function(x,y) {
  r <- x + (y-1)*(y-2)/2
  return(r)
}

mattovec <- function(mat){
  vec <- mat[upper.tri(mat, diag = FALSE)]
  return(vec)
}

log_normalization = function(x){

  x <- as.matrix(x)

  n <- dim(x)[2]

  gene.exprs.count <- rowSums(x != 0)

  sf <- colSums(x)/median(colSums(x))

  return(log(sweep(x, 2, sf, '/')+1))

}


keras_lasso_regression <- function(Y, X,  epochs = 100, batch_size = 128, lambda1 = 1.0, lambda2 = 1e10,
                                   learning_rate = 0.0001, lasso_threshold = 0, verbose = TRUE){

  p <- ncol(X)
  n <- nrow(X)

  X <- array_reshape(X, dim = c(nrow(X), ncol(X)))
  Y <- array_reshape(Y, dim = c(nrow(X), ncol(X)))


  network <- keras::keras_model_sequential() %>%
    layer_dense(units = p, input_shape = p, activation = "linear", kernel_regularizer = function(weight_matrix)
      regularizer_define(weight_matrix = weight_matrix ,lambda1 = lambda1,  lambda2 = lambda2),
      use_bias = FALSE, kernel_constraint = constraint_nonneg(), kernel_initializer = "zeros")


  network %>% keras::compile(loss = loss_define, optimizer = optimizer_rmsprop(lr = learning_rate), metrics = list("mean_squared_error"))
  history <- network %>% keras::fit(x = X, y = Y, epochs = epochs, batch_size = batch_size, validation_split = 0,
                                    verbose = verbose)


  weight <- keras::get_weights(network)[[1]]

  weight_final <- weight

  weight_final[weight_final <= lasso_threshold] <- 0

  return(weight_final)
}



regularizer_define <- function(weight_matrix ,lambda1 = 1.0, lambda2 = 1e10){
  lambda1 * k_sum(k_abs(weight_matrix), axis = c(1,2)) + lambda2 * tf$linalg$trace(tf$square(weight_matrix))
}


objective <- function(X, lambda, A, S){

  #objval <- (0.5*norm(X - A%*%X - X%*%S, 'F')^2 + lambda*sum(abs(A)) + lambda*sum(abs(S)))
  objval <- (0.5*norm(X - A%*%X - X%*%S, 'F')^2 + lambda*sum(abs(S)))
  return(objval)
}


loss_define <- function(y_true, y_pred){
  0.5 * k_sum(k_square(y_true - y_pred), axis = c(1,2))
}


calc.lambda <- function(X_count){

  X <- log_normalization(X_count)

  lambda <- sd(X)

}




calc.size.factor <- function(x) {

  sf <- Matrix::colSums(x)/median(Matrix::colSums(x))

  scale.sf <- 1

  list(unname(sf), scale.sf)

}




clean.data <- function(x) {

  if (!sum(grepl("matrix", class(x), ignore.case = TRUE))) {

    x <- Matrix::Matrix(as.matrix(x))

    message("Converting x to matrix.")

    if (!is.numeric(x)) {

      warning("Make sure x is numeric.")

    }

  }

  np <- dim(x)

  size <- as.numeric(np[1])*as.numeric(np[2])

  if(size > 2^31-1){

    inds <- split(1:np[2], ceiling(1:np[2]/1000))

    for(i in 1:length(inds)){

      x[, inds[[i]]][x[, inds[[i]]] < 0.001] <- 0

    }

  } else {

    x[x < 0.001] <- 0

  }

  if (is.null(np) | (np[2] <= 1))

    stop("x should be a matrix with 2 or more columns")

  if (min(Matrix::colSums(x)) == 0) {

    nzerocells <- sum(Matrix::colSums(x) == 0)

    x <- x[, Matrix::colSums(x) != 0]

    message("Removing ", nzerocells, " cell(s) with zero expression.")

  }

  if (is.null(rownames(x))) {

    rownames(x) <- 1:np[1]

  }

  x

}



fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright",
                          "left", "center", "right",
                          "bottomleft", "bottom", "bottomright"))

  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    }
  }

  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100

  x1 <- switch(pos,
               topleft     =x[1] + sw,
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)

  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)

  old.par <- par(xpd=NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}





