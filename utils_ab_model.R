#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Functions for the agent-based model
# Date: 30/07/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

## Functions ----
#--------------------------------------------------------------------------#

#' Compute pairwise distances between points in 2D or 3D space.
#'
#' This function computes pairwise Euclidean distances between points in a 2D or 3D space 
#' defined by the \code{x}, \code{y}, and optionally \code{z} coordinates of each point.
#'
#' @param dataframe A tibble or data frame containing at least the following columns:
#'   \code{x}: X coordinates of points,
#'   \code{y}: Y coordinates of points,
#'   \code{img_name}: Image name or identifier associated with each point.
#'   If \code{z} coordinates are provided, include a column named \code{z}.
#'   
#' @param use_z Logical, indicating whether to include the z dimension in distance calculation.
#'   Default is \code{FALSE}.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{img_name}{Image name associated with each pair of points.}
#'     \item{p1, p2}{Indices of the points used to compute the distance.}
#'     \item{dist}{Euclidean distance between the points specified by \code{p1} and \code{p2}.}
#'   }

compute_pairwise_distances <- function(df, z_dim = FALSE) {
  # Number of points
  n <- nrow(df)
  
  # Prepare result storage
  distances <- expand_grid(p1 = 1:n, p2 = 1:n) %>% # generate of possible pairs
    filter(p1 < p2) %>%  # avoid redundancy + self-pairing
    mutate(img_name = df$img_name[p1]) # add image name
  
  if (!z_dim) {
    # 2D case
    distances <- distances %>% 
      mutate( 
        # compute distance
        dist = sqrt(
          (df$x[p2] - df$x[p1])^2 + 
            (df$y[p2] - df$y[p1])^2
        )
      ) %>% 
      dplyr::select(img_name, p1, p2, dist)
    
  } else {
    # 3D case
    distances <- distances %>% 
      mutate( 
        # compute distance
        dist = sqrt(
          (df$x[p2] - df$x[p1])^2 + 
            (df$y[p2] - df$y[p1])^2 +
            (df$z[p2] - df$z[p1])^2
        )
      ) %>% 
      dplyr::select(img_name, p1, p2, dist)
  }
  
  return(distances)
}

#' Calculate 3D kernel density estimate with customizable bandwidth and grid size.
#'
#' This function computes a kernel density estimate (KDE) in three dimensions (3D)
#' using the provided data frame with x, y, and z coordinates.
#'
#' @param df A data frame containing the following columns:
#'   \describe{
#'     \item{x}{Numeric vector of x-coordinates.}
#'     \item{y}{Numeric vector of y-coordinates.}
#'     \item{z}{Numeric vector of z-coordinates.}
#'   }
#' @param h Scalar bandwidth for the KDE, which is used to create a diagonal 
#'   bandwidth matrix.
#' @param vol A named list containing the volume dimensions with components:
#'   \describe{
#'     \item{x}{Length of the volume in the x-direction.}
#'     \item{y}{Length of the volume in the y-direction.}
#'     \item{z}{Length of the volume in the z-direction.}
#'   }
#' @param d_res Resolution for the grid on which the KDE is evaluated. Default is 64.
#'   
#' @return A list with kernel density estimate results.
calculate_density_3d <- function(df, h, vol, d_res = 64) {
  
  # Matrix bandwidth from scalar bandwidth
  H <-  matrix(c(h, 0, 0, 0, h, 0, 0, 0, h), ncol = 3)
  
  # Output grid from image volume
  m <- crossing(
    x = seq(from = 0, to = vol$x, by = d_res),
    y = seq(from = 0, to = vol$y, by = d_res),
    z = seq(from = 0, to = vol$z, by = d_res)
  ) %>% 
    as.matrix()
  
  # Compute KDE
  kde <- ks::kde(x = as.matrix(df %>% dplyr::select(x, y, z)), H = H, eval.points = m)
  
  return(kde)
}


#' Calculate gradient of 3D kernel density estimate.
#'
#' This function computes the gradient of a 3D kernel density estimate (KDE) using
#' centered finite differences. It takes a KDE object as input and returns the gradient
#' components (dx, dy, dz) at each grid point.
#'
#' @param kde A list object returned by \code{\link[kedd]{kde}} function, containing:
#'   \describe{
#'     \item{eval.points}{A list of length 3 containing grid coordinates in each dimension (x, y, z).}
#'     \item{estimate}{A 3D array of density values corresponding to the grid coordinates.}
#'   }
#'   
#' @return A list with gradient components dx, dy, and dz, each being a 3D array of the same dimensions
#'   as the density estimate grid.
calculate_density_gradient_3d <- function(kde) {
  # Retrieve grid coordinates
  #x <- kde$eval.points[[1]]
  #y <- kde$eval.points[[2]]
  #z <- kde$eval.points[[3]]
  eval_points <- as_tibble(kde$eval.points)
  x <- eval_points %>% distinct(x) %>% pull(x)
  y <- eval_points %>% distinct(y) %>% pull(y)
  z <- eval_points %>% distinct(z) %>% pull(z)
  
  # Reformat density values as a 3D array
  d <- array(
    kde$estimate, 
    dim = c(
      length(z),
      length(y),
      length(x)
    ))
  # And reorder dimensions
  d <- aperm(d, c(3, 2, 1))
  
  #x <- unique(as.data.frame(kde$eval.points)$x)
  #y <- unique(as.data.frame(kde$eval.points)$y)
  #z <- unique(as.data.frame(kde$eval.points)$z)
  #d <- kde$estimate
  
  # Initialize arrays for the gradients with NA
  dx <- array(NA, dim = dim(d))
  dy <- array(NA, dim = dim(d))
  dz <- array(NA, dim = dim(d))
  
  # Calculate the gradient at each grid point using centered finite differences
  for (i in 2:(length(x) - 1)) {
    for (j in 2:(length(y) - 1)) {
      for (k in 2:(length(z) - 1)) {
        dx[i, j, k] <- (d[i + 1, j, k] - d[i - 1, j, k]) / (x[i + 1] - x[i - 1])
        dy[i, j, k] <- (d[i, j + 1, k] - d[i, j - 1, k]) / (y[j + 1] - y[j - 1])
        dz[i, j, k] <- (d[i, j, k + 1] - d[i, j, k - 1]) / (z[k + 1] - z[k - 1])
      }
    }
  }
  
  # Return the gradients as a list
  return(list(dx = dx, dy = dy, dz = dz))
}


#' Extract gradient values at given points from 3D kernel density estimate.
#'
#' This function extracts gradient values (dx, dy, dz) at specified points from a 3D
#' kernel density estimate (KDE) computed using the provided dataframe and gradient
#' components.
#'
#' @param df A data frame containing the following columns:
#'   \describe{
#'     \item{x}{Numeric vector of x-coordinates of points of interest.}
#'     \item{y}{Numeric vector of y-coordinates of points of interest.}
#'     \item{z}{Numeric vector of z-coordinates of points of interest.}
#'   }
#' @param gradient A list containing gradient components obtained from \code{calculate_density_gradient_3d}.
#'   Should include elements:
#'   \describe{
#'     \item{dx}{3D array of gradient values in the x-direction.}
#'     \item{dy}{3D array of gradient values in the y-direction.}
#'     \item{dz}{3D array of gradient values in the z-direction.}
#'   }
#' @param kde A list object returned by \code{calculate_density_3d} function, containing:
#'   \describe{
#'     \item{eval.points}{A list of length 3 containing grid coordinates in each dimension (x, y, z).}
#'   }
#'   
#' @return A list with gradient values (dx, dy, dz) at each specified point.
extract_gradient <- function(df, gradient, kde) {
  
  # Extract x, y, and z coordinates of original points
  x <- df$x
  y <- df$y
  z <- df$z
  
  # Extract x, y and z coordinates of density grid
  eval_points <- as_tibble(kde$eval.points)
  x_eval <- eval_points %>% distinct(x) %>% pull(x)
  y_eval <- eval_points %>% distinct(y) %>% pull(y)
  z_eval <- eval_points %>% distinct(z) %>% pull(z)
  
  # Initialize storage
  gradient_at_points <- list(
    dx = numeric(length(x)),
    dy = numeric(length(y)),
    dz = numeric(length(z))
  )
  
  # Find closest grid points to each input point
  for (i in 1:length(x)) {
    # Find indices of closest grid points in each dimension
    idx_x <- findInterval(x[i], x_eval)
    idx_y <- findInterval(y[i], y_eval)
    idx_z <- findInterval(z[i], z_eval)
    
    # Extract gradient values at closest grid points
    gradient_at_points$dx[i] <- gradient$dx[idx_x, idx_y, idx_z]
    gradient_at_points$dy[i] <- gradient$dy[idx_x, idx_y, idx_z]
    gradient_at_points$dz[i] <- gradient$dz[idx_x, idx_y, idx_z]
  }
  
  # Replace missing gradients by 0
  gradient_at_points$dx[is.na(gradient_at_points$dx)] <- 0
  gradient_at_points$dy[is.na(gradient_at_points$dy)] <- 0
  gradient_at_points$dz[is.na(gradient_at_points$dz)] <- 0
  
  return(gradient_at_points)
}


#' Update point positions based on density gradient and displacement length, for a specified proportion of points.
#'
#' This function updates the positions of points based on the gradient of a 3D
#' density estimate and a specified displacement length. It calculates the normalized
#' gradients (dx_norm, dy_norm, dz_norm) from the provided gradient values and then
#' randomly selects a proportion of points to update along these gradients.
#'
#' @param df A data frame containing the following columns:
#'   \describe{
#'     \item{x}{Numeric vector of x-coordinates of points to be updated.}
#'     \item{y}{Numeric vector of y-coordinates of points to be updated.}
#'     \item{z}{Numeric vector of z-coordinates of points to be updated.}
#'     \item{img_name}{Character vector specifying the image name associated with each point.}
#'   }
#' @param gradient_values A list containing gradient components obtained from \code{calculate_density_gradient_3d}.
#'   Should include elements:
#'   \describe{
#'     \item{dx}{Numeric vector or array of gradient values in the x-direction.}
#'     \item{dy}{Numeric vector or array of gradient values in the y-direction.}
#'     \item{dz}{Numeric vector or array of gradient values in the z-direction.}
#'   }
#' @param d_length Numeric scalar specifying the displacement length for updating points.
#' @param prop_mv Numeric scalar specifying the proportion of points to update. Should be between 0 and 1.
#'
#' @return A data frame with updated point positions (x_new, y_new, z_new) and corresponding image names (img_name).
move_points <- function(df, gradient_values, d_length, prop_mv = 1) {
  # Extract x, y, and z coordinates from dataframe
  x <- df$x
  y <- df$y
  z <- df$z
  
  # Extract gradient values
  dx <- gradient_values$dx
  dy <- gradient_values$dy
  dz <- gradient_values$dz
  
  # Calculate magnitudes of gradients
  magnitudes <- sqrt(dx^2 + dy^2 + dz^2)
  
  # Ensure magnitudes are not zero to avoid division by zero
  magnitudes[magnitudes == 0] <- 1
  
  # Normalize gradients
  dx_norm <- dx / magnitudes
  dy_norm <- dy / magnitudes
  dz_norm <- dz / magnitudes
  
  # Determine number of points to update based on proportion
  num_points <- nrow(df)
  num_to_update <- round(prop_mv * num_points)
  
  # Randomly select indices of points to update
  indices_to_update <- sample(num_points, size = num_to_update, replace = FALSE)
  
  # Update selected point positions based on density gradient and displacement length
  x_new <- x
  y_new <- y
  z_new <- z
  
  x_new[indices_to_update] <- x[indices_to_update] + d_length * dx_norm[indices_to_update]
  y_new[indices_to_update] <- y[indices_to_update] + d_length * dy_norm[indices_to_update]
  z_new[indices_to_update] <- z[indices_to_update] + d_length * dz_norm[indices_to_update]
  
  # Create a new dataframe with updated point positions and image names
  new_df <- tibble(img_name = df$img_name, x = x_new, y = y_new, z = z_new)
  
  return(new_df)
}



# Function to process data with given parameters
process_dataframe <- function(df, h = h, vol = vol, d_length, prop_mv) {
  # Compute distances before
  dist_before <- compute_pairwise_distances(df) %>% mutate(when = "before")
  
  # Compute 3D density
  dens_3d <- calculate_density_3d(df, h = h, vol = vol)
  
  # Compute density gradient
  grad_3d <- calculate_density_gradient_3d(dens_3d)
  
  # Extract gradient at points
  d_grad <- extract_gradient(
    df,
    gradient = grad_3d,
    kde = dens_3d
  )
  
  # Move points
  new_d_points <- move_points(
    df,
    gradient_values = d_grad,
    d_length = d_length,
    prop_mv = prop_mv
  )
  
  # Compute distances after moving
  dist_after <- compute_pairwise_distances(new_d_points) %>% mutate(when = "after")
  
  # Return distances before and after
  bind_rows(dist_before, dist_after)
}