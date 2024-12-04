#' Crop an image around a centerpoint
#' @export
crop_cimg_around_xy <- function(img, x, y, w=50){
  
  # Define cropping coordinates
  x_min <- x - w/2 |> round()
  x_max <- x + w/2 |> round()
  y_min <- y - w/2 |> round()
  y_max <- y + w/2 |> round()
  
  # Crop the image
  cropped_img <- img[x_min:x_max, y_min:y_max,, ] |> 
    as.cimg()
  
  return(cropped_img)
}

#' A function to shift zero-frequency component to the center of the spectrum.
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
#' @export
fftshift <- function(x) {
  # Shift the zero-frequency component to the center of the spectrum
  n <- length(x)
  return(c(x[(n/2 + 1):n], x[1:(n/2)]))
}

#' FFT Based Cross Correlation Computations
#' 
#' - http://www.learnpiv.org/fft/
#'   - https://chatgpt.com/share/674f5516-bb78-800f-9038-e00152429327
#'
#' Define a function for cross-correlation
#' @export
cross_correlation <- function(img1, img2) {
  # Ensure both images are matrices of the same size
  if (!all(dim(img1) == dim(img2))) {
    stop("Images must have the same dimensions for cross-correlation.")
  }
  
  # Compute cross-correlation using fft (fast Fourier transform)
  # Step 1: Compute FFT of both images
  fft1 <- fft(img1)
  fft2 <- fft(img2)
  
  # Step 2: Element-wise multiplication of one FFT with the conjugate of the other
  fft_cross <- fft1 * Conj(fft2)
  
  # Step 3: Inverse FFT to get cross-correlation
  cross_corr <- Re(fft(fft_cross, inverse = TRUE))
  
  # Step 4: Apply fftshift to center the zero-frequency component
  cross_corr_shift <- fftshift(cross_corr)
  
  # Convert back to a cimg object.
  cross_corr_shift_img <- cross_corr_shift |>
    array(dim(cross_corr)) |> 
    as.cimg()
  
  # Normalize the result
  cross_corr_norm <- cross_corr_shift_img / max(abs(cross_corr_shift))
  
  return(cross_corr_norm)
}

#' Compute image alignment offsets from FFT-based cross-correlation
#' @export
img_alignment_offset <- function(img1, img2){
  # Compute cross-correlation
  cross_corr <- cross_correlation(img1, img2)
  
  # Display the result
  # result |> as.cimg() |> plot()
  # image(cross_corr[,,1,], main = "Cross-Correlation", col = heat.colors(256))
  
  max_index <- which(cross_corr == max(cross_corr), arr.ind = TRUE) 
  
  # # Subtract the center of the image, because the correlation matrix may be double
  # # the size of the input image (padded to compute the correlation).
  # offset_y <- max_index[1] - nrow(img1) / 2
  # offset_x <- max_index[2] - ncol(img1) / 2
  # 
  # # Shift the second image to align with the first
  # aligned_image <- imager::imshift(img2, x = -offset_x, y = -offset_y)
  
  f <- function(o, d) if(o > d) {o - d} else {o + d}
  
  # Get pixel indexes of maximum cross-correlation.
  offsets <- max_index[1:2]
  # Remainder of division (deprecated after zero-centering in cross_correlation).
  # offsets <- max_index[1:2] %% dim(img1)[1:2] - 1
  
  offset_x <- offsets[1] - round(dim(img1)[1]/2)
  offset_y <- offsets[2] - round(dim(img1)[2]/2)
  cat("Offset X:", offset_x, "Offset Y:", offset_y, "\n")
  
  return(
    c(offset_x=offset_x,
      offset_y=offset_y,
      max_cor_x = offsets[1],
      max_cor_y = offsets[2],
      max_index=max_index[1:2]
      )
  )
}
