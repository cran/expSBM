#
#============= Initialize VEM algorithm via node allocation
#

expSBM_init <- function(edgelist, K, soft = TRUE)
{
  N <- max( c(edgelist[,1],edgelist[,2]) )
  nodes <- 1:N
  xint <- xnoint <- rep(0, N)
  matu <- matrix(0, N, N)
  
  int <- (edgelist[,3] == 1)
  noint <- (edgelist[,3] == 0)
  
  for ( i in 1:N ) {
    
    # start node allocation
    setu <- (edgelist[,1] == i)
    y <- edgelist[setu & int, c(2,4), drop = FALSE]
    w <- edgelist[setu & noint, c(2,4), drop = FALSE]
    iy <- unique(y[,1])
    iw <- unique(w[,1])
    xint[iy] <- tapply(y[,2], INDEX = y[,1], sum)
    # xint[iw] <- tapply(w[,2], INDEX = w[,1], mean)
    matu[i, nodes] <- xint
    
    xint[] <- 0
    xnoint[] <- 0
  }
  
  # lmatu <- ifelse(matu > 0, log(matu), 0)
  # lmatu <- ifelse( matu > 0, (matu - min(matu[matu!=0]))/(max(matu) - min(matu[matu!=0])), 0 )
  # z <- cmdscale(dist(lmatu))
  # z <- cmdscale( 1 - (lmatu + t(lmatu))/2 )
  
  matu <- (matu + t(matu))/2
  lmatu <- ifelse(matu > 0, log(matu), 0)
  A <- (lmatu + t(lmatu))/2
  d <- apply(A, 1, sum) + 1
  L <- diag(d) - A                            # unnormalized version
  L <- diag(d^-0.5) %*% L %*% diag(d^-0.5)    # normalized version
  ev <- eigen(L, symmetric = TRUE)
  # xx <- ev$vectors[,(ncol(ev$vectors)-K+1):ncol(ev$vectors)]    # find a lot of spurious solutions
  xx <- ev$vectors[,(ncol(ev$vectors)-5):ncol(ev$vectors)]
  fit <- kmeans(xx, centers = K, nstart = 500)
  zz <- mclust::unmap(fit$cluster) # + matrix(runif(N*K, 0,0.1), N,K)
  zz <- sweep(zz, 1, rowSums(zz), "/")
  classi <- mclust::unmap(fit$cluster)
  # 
  # if ( is.null(dim(xx)) ) { 
  #   fit <- Mclust(xx, G = K, prior = priorControl(), verbose = FALSE) 
  # } else { fit <- Mclust(xx, G = K, modelNames = c("EII", "VII", "EEI", "EVI", "VEI", "VVI"), prior = priorControl(), verbose = FALSE) }
  # zz <- fit$z + matrix(runif(N*K, 0,0.2), N,K)
  # zz <- sweep(zz, 1, rowSums(zz), "/")
  # classi <- mclust::unmap(fit$classification)
  
  return( list(Z = if ( soft ) zz else classi) )
}
