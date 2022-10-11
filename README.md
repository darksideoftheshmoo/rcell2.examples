# rcell2 examples

Examples for [rcell2](https://github.com/darksideoftheshmoo/rcell2) and companion packages, [rcell2.cellid](https://github.com/darksideoftheshmoo/rcell2-cellid) and [rcell2.magick](https://github.com/darksideoftheshmoo/rcell2-magick).

The `rcell2` packages provide functions to analyze Cell-ID single-cell cytometry data in a tidy and shiny framework.

It also includes sample images to test out the `rcell2` packages, located at `inst/extdata`.

# Installation

From R, run:

```r
if(!require(remotes)){
  # Install the "remotes" package if not found:
  install.packages("remotes")
}

# Install rcell2.examples
remotes::install_github("darksideoftheshmoo/rcell2.examples")
```

This package has no dependencies for installation.

However, all example Rmd notebooks use the rcell2 packages, and many other packages are required to run all the example notebooks.

It is left to the user to install them as needed.

# Development notebooks

Example data and analysis notebooks can be found in the examples package: [`rcell2.examples`](https://github.com/darksideoftheshmoo/rcell2.examples/tree/main).

The `inst/testings` directory holds many rmarkdown notebooks, where we explore different analysis approaches to single cell images and cytometry:

* Spatial distribution of fluorescent signals.
* Pattern detection in time-series.
* Classification examples.
* Cell boundary curvature analysis and alignment.
* ...

!['Automating' comes from the roots 'auto-' meaning 'self-', and 'mating', meaning 'screwing'.](https://imgs.xkcd.com/comics/automation.png)

# Todo

* Cleanup and test all notebooks.
