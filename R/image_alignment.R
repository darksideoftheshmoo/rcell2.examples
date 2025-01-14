#' Crop an image around a centerpoint
#' @export
crop_cimg_around_xy <- function(img, x, y, w=50){
  
  # Define cropping coordinates
  x_min <- x - w/2 |> round()
  x_max <- x + w/2 |> round()
  
  y_min <- y - w/2 |> round()
  y_max <- y + w/2 |> round()
  
  # Crop the image.
  # cropped_img <- img[y_min:y_max, x_min:x_max,, ] |> as.cimg()
  
  # Crop the image in the WRONG WAY.
  # Who knows why this makes it work with the other nasty hack.
  cropped_img <- img[x_min:x_max, y_min:y_max,, ] |> as.cimg()
  
  return(cropped_img)
}

#' Perform a Circular Shift on an Image
#'
#' Shifts the input image by half its dimensions along both axes, wrapping the image 
#' content around to ensure continuity at the boundaries. The function uses Fourier-style 
#' boundary conditions for wrapping.
#'
#' @param img A 2D or 3D array representing the image to be shifted. The dimensions 
#' should be specified as rows and columns for proper operation.
#'
#' @return A shifted image array where the content has been wrapped around to maintain 
#' continuity. The size and type of the returned array match the input.
#'
#' @details The function calculates the shift amounts as half of the image's dimensions 
#' (floored) along each axis. It then applies the shift using the `imshift` function with 
#' `boundary_conditions = 2`, which corresponds to a Fourier-style wrapping of the image 
#' content at the edges.
#'
#' @examples
#' # Create a sample matrix
#' img <- matrix(1:9, nrow = 3, ncol = 3)
#'
#' # Apply a circular shift
#' shifted_img <- circular_shift(img)
#'
#' @seealso \code{\link[imager]{imshift}} for more details on shifting and boundary conditions.
#'
#' @export
circular_shift <- function(img) {
  dims <- dim(img)
  shift_x <- floor(dims[2] / 2)
  shift_y <- floor(dims[1] / 2)
  shifted <- imshift(img, delta_x = shift_x, delta_y = shift_y, boundary_conditions = 2)
  return(shifted)
}

#' Shift the Zero-Frequency Component to the Center of the Spectrum
#'
#' Reorders the elements of a 1D or multi-dimensional array such that the zero-frequency 
#' component is moved to the center of the spectrum. This operation is commonly used in 
#' Fourier transform processing to improve visualization and analysis.
#'
#' @param x A numeric vector, matrix, or multi-dimensional array. The input represents 
#' a spectrum or image where the zero-frequency component is to be shifted.
#'
#' @return A numeric array of the same dimensions as the input, with the zero-frequency 
#' component shifted to the center. The data ordering is circularly wrapped.
#'
#' @details
#' The function divides the input into two halves and swaps them to place the zero-frequency 
#' component at the center. The shift is computed based on the input length and is effective 
#' for 1D data. For multi-dimensional data, the function can be extended or applied separately 
#' to each dimension.
#'
#' Compared to \code{\link{circular_shift}}, which operates on spatial domain images, 
#' \code{fftshift} focuses on rearranging frequency domain data for analysis. The key 
#' distinction is that \code{fftshift} is explicitly designed for frequency domain applications, 
#' whereas \code{circular_shift} applies a spatial circular wrap using Fourier-style boundary 
#' conditions.
#' 
#' - http://matlab.izmiran.ru/help/techdoc/ref/fftshift.html
#' - https://stackoverflow.com/a/49863001
#' 
#' >  FFTSHIFT shifts the zero-frequency component to the center of the signal.
#' > In this case the signal is an image. A good visual guide is this.
#' > If you expand the original output image, you will see something akin to this:
#'
#' Esto fue importante para poder encontrar consistentemente el máximo de la 
#' correlacion cruzada cerca del centro si el offset es pequeño.
#'
#' @examples
#' # Shift a 1D array
#' vec <- 1:8
#' fftshift(vec)
#'
#' # Shift a 2D matrix (row and column shift must be implemented separately)
#' mat <- matrix(1:16, nrow = 4, ncol = 4)
#' row_shifted <- fftshift(mat[1, ])
#'
#' @seealso \code{\link{circular_shift}} for spatial domain circular shifting.
#'
#' @export
fftshift <- function(x) {
  # Shift the zero-frequency component to the center of the spectrum
  n <- length(x)
  x_shifted <- c(x[(n/2 + 1):n], x[1:(n/2)])
  result <- x_shifted |> array(dim(x))
  return(result)
}

#' FFT Based Cross Correlation Computations
#' 
#' - http://www.learnpiv.org/fft/
#'   - https://chatgpt.com/share/674f5516-bb78-800f-9038-e00152429327
#'
#' Define a function for cross-correlation
#' @export
cross_correlation <- function(img1, img2, zero_shift=F, circ_shift=T) {
  # Ensure both images are matrices of the same size
  if (!all(dim(img1) == dim(img2))) {
    stop("Images must have the same dimensions for cross-correlation.")
  }
  
  # Compute cross-correlation using fft (fast Fourier transform).
  # Step 1: Compute FFT of both images.
  fft1 <- fft(img1)
  fft2 <- fft(img2)
  
  # Step 2: Compute Cross Power Spectrum.
  # Element-wise multiplication of one FFT with the conjugate of the other.
  fft_cross <- fft1 * Conj(fft2)
  
  # Step 3: Inverse FFT to get cross-correlation.
  phase_corr <- Re(fft(fft_cross, inverse = TRUE))

  # Normalize by the number of pixels
  phase_corr <- Re(phase_corr) / prod(dim(img1))
  
  # Apply FFT-shift to center the zero-frequency component.
  if(zero_shift) fft_cross <- fftshift(phase_corr) |> as.cimg()
  
  # Apply circular-shift to center the zero-frequency component.
  if(circ_shift) phase_corr <- circular_shift(phase_corr)
  
  # Normalize the result.
  cross_corr_norm <- phase_corr / max(abs(phase_corr))
  
  return(cross_corr_norm)
}

get_value <- function(col_i, row_j, m) {
  dims <- dim(m)
  rows_j <- dims[1] # y = rows = j
  cols_i <- dims[2] # x = cols = i
  if (row_j > 0 && row_j <= rows_j && col_i > 0 && col_i <= cols_i) {
    return(m[row_j, col_i])
  } else {
    return(NA) # Boundary handling (zero/NA padding)
  }
}

#' Compute image alignment offsets from FFT-based cross-correlation
#' @export
img_alignment_offset <- function(img1, img2, s=NULL, s2=NULL, verbose=F){
  # Compute cross-correlation
  cross_corr <- cross_correlation(img1, img2)
  
  # Weigh before finding the max.
  if(!is_null(s)){
    if(verbose) cat("Applying gaussian weighing.")
    cross_corr[,,1,1] <- cross_corr[,,1,1] |> apply_gaussian_weighting(s=s)
  }
  if(!is.null(s2)){
    if(verbose) cat("Applying blur.")
    cross_corr <- cross_corr |> imager::isoblur(sigma=s2)
  }

  # Step 1: Find the peak in the phase correlation
  max_index <- which(cross_corr == max(cross_corr), arr.ind = TRUE)
  offsets <- max_index[1:2]
  
  # Get pixel indexes of maximum cross-correlation.
  # Subtracting one because I observe drift otherwise.
  # offset_y <- offsets[1] - round(dim(img1)[1]/2) # y
  # offset_x <- offsets[2] - round(dim(img1)[2]/2) # x
  
  # Step 2: Extract a 3x3 neighborhood around the peak
  peak_y <- offsets[1] # Row index.
  peak_x <- offsets[2] # Column index.
  neigh_n <- 2  # Amount of neighbours.
  neighborhood <- expand.grid(
    x = (peak_x - neigh_n):(peak_x + neigh_n),
    y = (peak_y - neigh_n):(peak_y + neigh_n)
  )
  neighborhood$z <- apply(neighborhood, 1, function(idx) {
    get_value(col_i = idx["x"],
              row_j = idx["y"], 
              m = cross_corr[,,1,1]) # Note: y comes first in R matrices
  })
  
  # Step 3: Fit a 2D parabola using lm.
  fit <- lm(
    z ~ poly(x, 2, raw = TRUE) + poly(y, 2, raw = TRUE) + x:y,
    data = neighborhood
  )
  
  # Step 4: Generate an over-sampled grid of sub-pixel coordinates.
  subpixel_range <- seq(-neigh_n, neigh_n, by=0.05)
  oversampled_grid <- expand.grid(
    x = peak_x + subpixel_range,
    y = peak_y + subpixel_range
  )
  
  # Step 5: Predict values on the over-sampled grid.
  oversampled_grid$z <- predict(fit, newdata = oversampled_grid)
  
  # Step 6: Find the maximum predicted value.
  center_x <- (dim(img1)[2] + 1) / 2
  center_y <- (dim(img1)[1] + 1) / 2
  max_pred_idx <- which.max(oversampled_grid$z)
  offset_y <- oversampled_grid$y[max_pred_idx] - center_y # y
  offset_x <- oversampled_grid$x[max_pred_idx] - center_x # x
  
  if(verbose) cat("Offset X:", offset_x, "Offset Y:", offset_y, "\n")
  
  return(
    c(offset_x=offset_x,
      offset_y=offset_y,
      max_cor_x = offsets[1],
      max_cor_y = offsets[2],
      max_index=max_index[1:2]
      )
  )
}

#' generate_gaussian_matrix
#' @export
generate_gaussian_matrix <- function(n_row, n_col, s) {
  # Create a grid of (x, y) coordinates centered at (0, 0)
  x <- seq(-floor(n_col / 2), floor((n_col - 1) / 2), length.out = n_col)
  y <- seq(-floor(n_row / 2), floor((n_row - 1) / 2), length.out = n_row)
  
  grid_x <- matrix(rep(x, n_row), nrow = n_row)
  grid_y <- matrix(rep(y, each = n_col), nrow = n_row)
  
  # Compute the Gaussian density
  gaussian_matrix <- exp(-(grid_x^2 + grid_y^2) / (2 * s^2))
  
  # Normalize so the maximum value is 1
  gaussian_matrix <- gaussian_matrix / max(gaussian_matrix)
  
  return(gaussian_matrix)
}

#' apply_gaussian_weighting
#' @export
apply_gaussian_weighting <- function(m, s) {
  # Get matrix dimensions
  n_row <- dim(m)[1]
  n_col <- dim(m)[2]
  
  weight_m <- generate_gaussian_matrix(n_row = n_row, n_col = n_col, s = s)
  
  # Apply the Gaussian weighting
  weighted_matrix <- m * weight_m
  
  return(weighted_matrix)
}


# Function to add a border to an image
add_border <- function(img, border_width = 1, border_color = 0) {
  img_height <- dim(img)[1]
  img_width <- dim(img)[2]
  
  # Create a new image with an additional border
  img_with_border <- matrix(border_color, 
                            nrow = img_height + 2 * border_width, 
                            ncol = img_width + 2 * border_width) |> 
    as.cimg()
  
  # Place the original image in the center of the new image
  img_with_border[(border_width + 1):(border_width + img_height),
                  (border_width + 1):(border_width + img_width),1,1] <- img
  
  return(img_with_border)
}

#' Create a 2D Tile from a List of Images
#'
#' This function arranges a list of images into a 2D grid (tile) with a specified number of rows, columns, 
#' or aspect ratio. It can handle incomplete rows by padding with empty (transparent) images.
#'
#' @param images A list of images, where each image is an object of class `cimg` (e.g., from the `imager` package).
#' @param n_cols Integer specifying the number of columns in the tile. If `NULL`, it is calculated automatically 
#'   based on `n_rows` or the total number of images and the `aspect_ratio`.
#' @param n_rows Integer specifying the number of rows in the tile. If `NULL`, it is calculated automatically 
#'   based on `n_cols` or the total number of images and the `aspect_ratio`.
#' @param aspect_ratio Numeric value indicating the desired width-to-height ratio of the grid. Only used 
#'   if both `n_cols` and `n_rows` are `NULL`.
#' @param verbose Logical, if `TRUE`, prints the number of rows and columns being used in the grid.
#'
#' @return A `cimg` object representing the 2D tiled image.
#'
#' @details If both `n_cols` and `n_rows` are provided, they take precedence, and the `aspect_ratio` is ignored. 
#' The function pads incomplete rows with empty images to maintain a uniform grid structure.
#'
#' @seealso [imager::as.cimg], [imager::plot]
#'
#' @import imager
#' @export
create_image_tile <- function(images, n_cols = NULL, n_rows = NULL, aspect_ratio = 1, verbose = FALSE) {
  # Get the dimensions of the images
  img_width <- dim(images[[1]])[2]
  img_height <- dim(images[[1]])[1]
  n_images <- length(images)
  
  # Calculate number of rows and columns
  if (is.null(n_cols) & is.null(n_rows)) {
    # If neither rows nor columns are specified, approximate a square grid
    n_rows <- ceiling(sqrt(n_images * aspect_ratio))
    n_cols <- ceiling(n_images / n_rows)
  } else if (is.null(n_cols)) {
    n_cols <- ceiling(n_images / n_rows)
  } else if (is.null(n_rows)) {
    n_rows <- ceiling(n_images / n_cols)
  }
  
  if (verbose) {
    cat("Arranging into a grid with", n_cols, "rows and", n_rows, "columns.\n")
  }
  
  # Initialize an empty list to hold the rows of the tiled image
  rows_list <- list()
  
  for (i in 1:n_cols) {
    # Get the images for the current row
    row_images <- images[((i - 1) * n_rows + 1):min(i * n_rows, n_images)]
    
    # If there are fewer images than columns in the last row, pad with empty images
    if (length(row_images) < n_rows) {
      padding <- lapply(1:(n_rows - length(row_images)), function(x) imager::as.cimg(matrix(NA, nrow = img_height, ncol = img_width)))
      row_images <- c(row_images, padding)
    }
    
    # Concatenate the row of images horizontally
    row_images <- row_images |> 
      lapply(add_border) |> 
      lapply(function(i) i[,,1,1])
    row_img <- Reduce(function(x, y) cbind(x, y), row_images) 
    
    # Append to the rows list
    rows_list[[i]] <- row_img |> as.cimg()
  }
  
  # Now concatenate all rows vertically
  rows_list <- rows_list |> lapply(function(i) i[,,1,1])
  tile_image <- Reduce(function(x, y) rbind(x, y), rows_list) |> 
    as.cimg()
  
  return(tile_image)
}

#' adjust_offset
#' @export
adjust_offset <- function(img_paths, xy, t.subject, t.query, w=100, s=250, s2=0.0, verbose = F){
  
  # Get image paths.
  path1 <- img_paths |> dplyr::filter(t.frame == t.subject) |> with(file)
  path2 <- img_paths |> dplyr::filter(t.frame == t.query) |> with(file)
  
  # Load two images
  frame1 <- load.image(path1)
  frame2 <- load.image(path2)
  
  # Crop frame1.
  img1 <- frame1 |> crop_cimg_around_xy(x = xy[1], y = xy[2], w = w)
  # Crop frame2.
  img2 <- frame2 |> crop_cimg_around_xy(x = xy[1], y = xy[2], w = w)
  
  # Calculate the offset through cross-correlation.
  offsets <- img_alignment_offset(img1, img2, s=s, s2=s2, verbose = verbose)
  
  # Get the final XY offsets.
  xy_offsets <- offsets[c("offset_x", "offset_y")]
  return(xy_offsets)
}

#' Propagate Offsets Across Image Frames
#'
#' This function computes and propagates XY offsets across a sequence of image frames, aligning them based on a 
#' specified reference frame and an initial alignment point.
#'
#' @param img_paths A data frame containing information about the image frames. Must include columns `t.frame` 
#'   for the time points of the frames. Additional columns `x_offset`, `y_offset`, and `t.ref` will be updated or added.
#' @param t.bud Numeric, the initial reference time frame to align all subsequent frames to.
#' @param xy Numeric vector of length 2, the initial XY alignment point.
#' @param w Numeric, the width of the region used for alignment in the image. Default is 50.
#' @param s Optional parameter passed to `adjust_offset` to control the alignment strategy (default is `NULL`).
#' @param s2 Optional parameter passed to `adjust_offset` for additional alignment control (default is `NULL`).
#' @param xy_adj Logical, if `TRUE`, offsets are adjusted relative to the accumulated offset of previous frames. 
#'   If `FALSE`, all frames are aligned relative to the initial alignment point. Default is `FALSE`.
#' @param verbose Logical, if `TRUE`, prints detailed progress and alignment information. Default is `FALSE`.
#'
#' @return A data frame with the same structure as `img_paths`, including updated columns for `x_offset`, `y_offset`, 
#'   and `t.ref`, which store the computed offsets and reference time frame for each image.
#'
#' @details
#' - If `xy_adj` is `FALSE`, all frames are aligned relative to the event frame (the `t.bud` time frame) using the 
#'   provided `xy` alignment point.
#' - If `xy_adj` is `TRUE`, each frame is aligned relative to the accumulated offsets of previous frames. 
#'   Note that this mode can sometimes cause instability or drift.
#' - The function relies on `adjust_offset`, which performs the actual alignment between frames.
#'
#' @seealso [adjust_offset]
#' @export
propagate_offsets <- function(img_paths, t.bud, xy, w=50, s=NULL, s2=NULL, xy_adj=F, verbose=F){
  # Check
  stopifnot(length(xy) == 2, paste("The 'xy' argument is not a vector of length 2:", xy)) 
  # Initialize.
  img_paths[1, c("x_offset", "y_offset", "t.ref")] <- c(0, 0, t.bud)
  # Early return.
  if(nrow(img_paths) < 2) return(img_paths)
  # Process.
  for(i in 2:nrow(img_paths)){
    
    # Choose center point for alignment.
    if(!xy_adj){
      # Align all to the event frame.
      xy_next = xy
    } else {
      # Align centered over the accumulated offset.
      # Doing this for some reason messes up stuff.
      xy_i <- img_paths[1:(i-1), c("x_offset", "y_offset")] |> colSums(na.rm = T)
      xy_next = xy + xy_i
    }

    # Align with the next frame.
    t.subject <- img_paths$t.frame[i-1]
    t.query <- img_paths$t.frame[i]

    # Align.
    xy_offsets <- adjust_offset(
      img_paths = img_paths, 
      xy = xy_next,
      t.subject = t.subject, t.query = t.query,
      w=w, s=s, s2=s2,
      verbose=verbose
    )
    
    # Append.
    img_paths[i, c("x_offset", "y_offset")] <- xy_offsets
    img_paths[i, "t.ref"] <- t.subject
  }
  
  return(img_paths)
}
