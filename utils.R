## Packages ----
#--------------------------------------------------------------------------#
suppressWarnings(library(tidyverse))
suppressWarnings(library(here))
suppressWarnings(library(reshape2))
suppressWarnings(library(parallel))

# Null hypothesis
suppressWarnings(library(MASS))
suppressWarnings(library(twosamples))
suppressWarnings(library(extraDistr))

# Reading data
suppressWarnings(library(arrow))

# Spatial point pattern
suppressWarnings(library(spatstat))
suppressWarnings(library(dixon))

# Co-occurrence network
suppressWarnings(library(cooccur))
suppressWarnings(library(visNetwork))

# Species association network
suppressWarnings(library(EMtree))
suppressWarnings(library(ade4))
suppressWarnings(library(PLNmodels))

# Plots
suppressWarnings(library(gridExtra))
suppressWarnings(library(ggraph))
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(tidygraph))
suppressWarnings(library(ggpubr))
suppressWarnings(library(patchwork))
suppressWarnings(library(cmocean))
suppressWarnings(library(scales))
suppressWarnings(library(ggtext))


theme_set(theme_minimal())

## Subsampling ----
#--------------------------------------------------------------------------#
sub_sample <- TRUE # whether to subsample or not
n_img <- 10000 # number of images to consider for subsampling
n_dist <- 10000 # number of distances to retain for comparisons


## Image volume ----
#--------------------------------------------------------------------------#
# image volume in pixels
vol <- c()
vol$x <- 10240
vol$y <- 2048
vol$z <- 9572

## Directories ----
#--------------------------------------------------------------------------#
data_dir <- here("data")

## Parallel ----
#--------------------------------------------------------------------------#
n_cores <- 12

## Seed ----
#--------------------------------------------------------------------------#
set.seed(1)

## Coast data ----
#--------------------------------------------------------------------------#
coast <- read_csv(file.path(data_dir, "raw/coast.csv"), show_col_types = FALSE)


## Functions ----
#--------------------------------------------------------------------------#

#' Compute distances between individuals of two taxonomic groups within an ensemble of images
#' 
#' In each image, the distance between all unique pairs of points is computed.
#'
#' @param points dataframe with columns `x` and `y` for point position, `img_name` for image reference and `taxon` for taxonomic groups
#' @param n_cores number of cores to use for parallel processing
#' @param z_dim whether to consider 3rd dimension (`z`)
#'
#' @return a dataframe with all distances for each image
#' 
compute_all_dist_inter <- function(points, n_cores = 12, z_dim = FALSE) {
  dist_all <- mclapply(unique(points$img_name), function(name) {
    
    # Get points within image, split dataframe by taxon
    if (z_dim){ # 3D case
      points_img <-  points %>% 
        filter(img_name == name) %>% 
        select(taxon, x, y, z) %>% 
        group_by(taxon) %>% 
        group_split(.keep = F)
    } else { # 2D case
      points_img <- points %>% 
        filter(img_name == name) %>% 
        select(taxon, x, y) %>% 
        group_by(taxon) %>% 
        group_split(.keep = F)
    }
    
    # Use proxy::dist to compute distances between two matrices
      dist_all <- proxy::dist(points_img[[1]], points_img[[2]])
      dist_all <- `dim<-`(c(dist_all), dim(dist_all))
      dist_all <- dist_all %>% 
        melt(varnames = c("p1", "p2")) %>% 
        as_tibble() %>% 
        rename(dist = value) %>% 
        mutate(img_name = name)
    
    return(dist_all)
   }, mc.cores = n_cores)
  
  # Combine results in a dataframe
  dist_all <- do.call(bind_rows, dist_all)

  # Check that we have the number of expected distances within each image
  # Number of distances in each image
  dist_ok <- dist_all %>% 
    count(img_name, name = "n_dist") %>% 
    # Join with theoretical number of distances for each image
    left_join(
      points %>% 
        count(img_name, taxon) %>% 
        group_by(img_name) %>% 
        summarise(n_dist_th = prod(n)),
      by = join_by(img_name)
    ) %>% 
    mutate(ok = n_dist == n_dist_th)
    
  if(!all(dist_ok$ok)){stop("Number of computed distances differs from what is expected based on the number of points and their type.")} # should return TRUE
  
  return(dist_all)
}




#' Compute distances between all points within an ensemble of images
#' 
#' In each image, the distance between all unique pairs of points is computed.
#'
#' @param points dataframe with columns `x` and `y` for point position and `img_name` for image reference
#' @param n_cores number of cores to use for parallel processing
#' @param z_dim whether to consider 3rd dimension (`z`)
#'
#' @return a dataframe with all distances for each image
#' 
compute_all_dist <- function(points, n_cores = 12, z_dim = FALSE) {
  dist_all <- mclapply(unique(points$img_name), function(name) {
    
    # Get points within image
    if (z_dim){ # 3D case
      points_img <-  points %>% 
        filter(img_name == name) %>% 
        select(x, y, z) %>% 
        as.matrix()
    } else { # 2D case
      points_img <-  points %>% 
        filter(img_name == name) %>% 
        select(x, y) %>% 
        as.matrix()
    }
    
    # Compute distances between points
    melt(as.matrix(dist(points_img)), varnames = c("p1", "p2")) %>% 
      as_tibble() %>% 
      filter(p1 != p2) %>% 
      rename(dist = value) %>% 
      mutate(img_name = name)
  }, mc.cores = n_cores)
  
  # Combine results in a dataframe
  dist_all <- do.call(bind_rows, dist_all)  %>% 
    # Keep only one of distances compute between A and B (A to B and B to A were computed)
    mutate(pair = ifelse(p1 < p2, paste(p1, p2), paste(p2, p1))) %>% 
    distinct(img_name, pair, dist, .keep_all = TRUE) %>% 
    select(-pair)
  
  # Check that we have the number of expected distances within each image
  # For a set of n points, the number of unique distances is n(n-1)/2
  dist_ok <- left_join(
    points %>% count(img_name, name = "n_obj"),
    dist_all %>% count(img_name, name = "n_dist"),
    by = join_by(img_name)
  ) %>% 
    mutate(ok = n_dist == (n_obj * (n_obj-1))/2)
  
  if(!all(dist_ok$ok)){stop("Number of computed distances differs from what is expected based on the number of points.")} # should return TRUE
  
  return(dist_all)
}

#' Resampling procedure for  edges probability
#'
#' @param counts Data of observed counts with dimensions n x p, either a matrix, data.frame or tibble.
#' @param covar_matrix matrix of covariates, should have the same number of rows as the count matrix.
#' @param unlinked An optional vector of nodes which are not linked with each other
#' @param O Matrix of offsets, with dimension n x p
#' @param user_covariance_estimation A user-provided function for the estimation of a covariance
#' @param v The proportion of observed data to be taken in each sub-sample. It is the ratio (sub-sample size)/n
#' @param S Total number of wanted sub-samples.
#' @param maxIter Maximum number of EMtree iterations at each sub-sampling.
#' @param cond.tol Tolerance for the psi matrix.
#' @param eps Precision parameter controlling the convergence of weights beta
#' @param cores Number of cores, can be greater than 1 if data involves less than about 32 species.
#' @param init boolean: should the resampling be carried out with different initial points (TRUE), or with different initial data (FALSE)
#' @return Returns a list which contains the Pmat data.frame, and vectors of EMtree maximum iterations and running times in each
#' resampling.
#'\itemize{
#'  \item{Pmat: }{S x p(p-1)/2 matrix with edge probabilities for each resample}
#'  \item{maxIter: }{EMtree maximum iterations in each resampling.}
#'  \item{times: }{EMtree running times in each resampling.}
#' }
#'
#' @export
#' @importFrom PLNmodels PLN
#' @importFrom parallel mclapply
#' @examples
#'n=100
#'p=12
#'S=5
#'set.seed(2021)
#'simu=data_from_scratch("erdos",p=p,n=n)
#'G=1*(simu$omega!=0) ; diag(G) = 0
#'# With default evaluation, using the PLNmodel paradigm:
#'default_resample=ResampleEMtree(simu$data, S=S,cores = 1)
#'
#'# With provided correlation estimation function:
#'estimSigma<-function(counts, covar_matrix, sample){
#'Dum_Sigma = cov2cor(cov(counts[sample,]))
#'}
#'custom_resample=ResampleEMtree(simu$data,S=S,cores = 1,user_covariance_estimation=estimSigma)
#'
#'# We then run the stability selection to find the optimal selection frequencies,
#'# for a stability of 85%:
#'stab_default=StATS(default_resample$Pmat, nlambda=50, stab.thresh=0.8,plot=TRUE)
#'stab_custom=StATS(custom_resample$Pmat, nlambda=50, stab.thresh=0.8,plot=TRUE)
#'
#' #Check quality of result
#'table(pred=1*(stab_default$freqs_opt>0.9), truth=ToVec(G))
#'table(pred=1*(stab_custom$freqs_opt>0.9), truth=ToVec(G))
ResampleEMtree <- function(counts,covar_matrix=NULL, unlinked=NULL,
                           O=NULL,user_covariance_estimation=NULL,
                           v=0.8, S=1e2, maxIter=30, cond.tol=1e-10,
                           eps=1e-3,cores=3, init=FALSE){
  cat("Computing",S,"probability matrices with", cores, "core(s)... ")
  t1<-Sys.time()
  counts<-as.matrix(counts)
  n <- nrow(counts);  p <- ncol(counts)
  P <- p * (p - 1) / 2 ; V <- round(v * n)
  Pmat <- matrix(0, S, P)
  
  #- offsets and covariates
  if(is.null(O)){ O<-matrix(1, n, p)}
  if(is.null(covar_matrix)){#default intercept
    X<-matrix(1,nrow=n,ncol=1)
  }else{X<-as.matrix(covar_matrix)}
  #- parallel computation of S fits of new_EMtree
  if(is.null(user_covariance_estimation)){
    suppressWarnings(
      PLNfit <- PLNmodels::PLN(counts ~ -1  + offset(log(O)) + .,
                               data=data.frame(X),control=PLN_param("trace"=0))
    )}
  obj<-parallel::mclapply(1:S,function(b){
    if(init){
      inf<-EMtree( PLNfit,unlinked,n=n, maxIter=maxIter, cond.tol=cond.tol,
                   verbatim=TRUE,eps=eps,plot=FALSE,
                   random.init = TRUE)[c("edges_prob","maxIter","timeEM")]
    }else{
      sample <- sample(1:n, V, replace = FALSE)
      if(is.null(user_covariance_estimation)){
        #inception
        M.sample<-PLNfit$var_par$M[sample,]
        S2.sample<-PLNfit$var_par$S2[sample,]
        CorY<-cov2cor(t(M.sample)%*%M.sample+diag(colSums(S2.sample)))
      }else{
        CorY<-user_covariance_estimation(counts=counts, covar_matrix=X,
                                         sample=sample)
      }
      try({
        inf<-EMtree( CorY,unlinked,n=n, maxIter=maxIter, cond.tol=cond.tol,
                     verbatim=TRUE,eps=eps,plot=FALSE,
                     random.init=TRUE)[c("edges_prob","maxIter","timeEM")]
      }, silent=TRUE)
      if(!exists("inf")) inf<-NA #depending on the sample drawn, it is possible that computation fail
      # because of bad conditioning of the Laplacian matrix of the weights beta.
      # This can happen especially when using the "unlinked" parameter.
    }
    return(inf)
  }, mc.cores=cores)
  bad_samples<-which(do.call(rbind, lapply(obj, length))!=3)
  time<-difftime(Sys.time(), t1)
  cat(round(time,2),  attr(time, "units"),"\n")
  if(length(bad_samples)!=0){
    cat(length(bad_samples), " failed samples.\n")
    obj<-obj[-bad_samples]
  }
  Pmat<-do.call(rbind,lapply(obj,function(x){ToVec(x$edges_prob)}))
  summaryiter <- do.call(c,lapply(obj,function(x){x$maxIter}))
  times<-do.call(c,lapply(obj,function(x){x$timeEM}))
  return(list(Pmat=Pmat,maxIter=summaryiter,times=times))
}