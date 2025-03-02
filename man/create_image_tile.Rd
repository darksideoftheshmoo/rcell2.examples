% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_alignment.R
\name{create_image_tile}
\alias{create_image_tile}
\title{Create a 2D Tile from a List of Images}
\usage{
create_image_tile(
  images,
  n_cols = NULL,
  n_rows = NULL,
  aspect_ratio = 1,
  verbose = FALSE
)
}
\arguments{
\item{images}{A list of images, where each image is an object of class `cimg` (e.g., from the `imager` package).}

\item{n_cols}{Integer specifying the number of columns in the tile. If `NULL`, it is calculated automatically 
based on `n_rows` or the total number of images and the `aspect_ratio`.}

\item{n_rows}{Integer specifying the number of rows in the tile. If `NULL`, it is calculated automatically 
based on `n_cols` or the total number of images and the `aspect_ratio`.}

\item{aspect_ratio}{Numeric value indicating the desired width-to-height ratio of the grid. Only used 
if both `n_cols` and `n_rows` are `NULL`.}

\item{verbose}{Logical, if `TRUE`, prints the number of rows and columns being used in the grid.}
}
\value{
A `cimg` object representing the 2D tiled image.
}
\description{
This function arranges a list of images into a 2D grid (tile) with a specified number of rows, columns, 
or aspect ratio. It can handle incomplete rows by padding with empty (transparent) images.
}
\details{
If both `n_cols` and `n_rows` are provided, they take precedence, and the `aspect_ratio` is ignored. 
The function pads incomplete rows with empty images to maintain a uniform grid structure.
}
\seealso{
[imager::as.cimg], [imager::plot]
}
