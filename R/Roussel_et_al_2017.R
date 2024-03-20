#' The recursive function H_n^k
#'
#' This function is the recursive function described in equation 7 (see references). The user does
#' not have to use this function. Function M provides the pre-processing procedure to compute H.
#' Note that R is very inefficient with recursive functions. This function should be written in C or
#' C++ to be efficient. However, for compatibility safety this function is written directly in R.
#'
#' @param n scalar. The bin number
#' @param P vector. The probabilities P_k^n without the null probabilities (see equation 4)
#' @param h vector. The heights h associated to P_k^n  probabilities (see equation 7)
#'
#' @references Roussel, J.-R., Caspersen, J., Béland, M., Thomas, S., & Achim, A. (2017). Removing
#' bias from LiDAR-based estimates of canopy height: Accounting for the effects of pulse density and
#' footprint size. Remote Sensing of Environment, 198, 1–16. https://doi.org/10.1016/j.rse.2017.05.032
#' @noRd
NULL

# ============= MOVED TO C++ ===============
#H = function(n, P, h)
#{
#  if(n > 1)
#    return( P[n]*h[n]+(1-P[n])*H(n-1, P, h) )
#  else
#    return( h[1] )
#}

#' Pre-processing function
#'
#' The pre-processing function which filters the null p_i probabilities and computes
#' the P_k^n probabilities described in equation 4.
#'
#' @param n scalar. The number of points
#' @param p vector. The probabilities p
#' @param h vector. The heights h associated to the probabilities
#'
#' @noRd
M = function(n, p, h)
{
  h = h[p>0]
  p = p[p>0]

  nbin = length(p)
  q = numeric(nbin)
  csum = cumsum(p)
  for(i in 1:nbin) q[i] = 1-(1-p[i]/csum[i])^n
  Hcpp(nbin, q, h)
}

#' Canopy histogram
#'
#' The standardized histogram of canopy (see figure 6 and 13) for one raster. This
#' function is applied in the function rasterize_metrics() of the lidR package.
#'
#' @param z vector. Z coordinates
#' @param zmax scalar. To be sure to capture the whole pattern set a value higher
#' than the highest value of the training dataset
#' @param bin scalar. The bin height to compute the histogram of probabilities
#'
#' @noRd
canopy_histogram = function(z, zmax, bin)
{
  z = z[z > 0]
  z = z-max(z)
  z = z[-which.max(z)]
  bk = seq(-zmax, 0, bin)
  hi = graphics::hist(z, breaks = bk, plot = F)

  h = hi$breaks[-1]
  p = as.list(hi$counts)
  names(p) = round(h,2)

  return(p)
}

#' Max height error estimation map
#'
#' This function evaluates, for a given pixel size the expected error made on the metrics hmax as a
#' function of the pulse density according to Roussel et al (2007) (see references). The metric hmax
#' is always underestimated because discrete laser pulses are statistically unlikely to hit the very
#' top of the canopy. The highest laser hit is necessarily below the theoretically highest point. How much?
#' It depends on the shape/complexity of the canopy and the sampling density. By evaluating the shape
#' of the canopy and using probability theory, it is possible to estimate the expected error made as
#' a function of the pulse density (see references).
#'
#' While the original paper estimate the "shape" of the canopy for a given point cloud and consequently
#' assigns an error for the overall point cloud, this functions estimates the "shape" of the canopy locally
#' using pixels larger that the requested resolution. For example if user want to know the error on hmax
#' for 20x20 m pixels, using a factor of 5, then, the function will estimate the error using a
#' 100x100 windows, i.e. 25 sub-pixel are used to estimate the error into the major 100x100 m pixel.
#'
#' @return A SpatRaster with a resolution f*res that tells, for each pixel the expected error on hmax
#' if hmax is computed in pixels of resolution 'res'.
#'
#' @param las LAS or LAScatalog from the lidR package. The point cloud must be normalized and
#' lake or any other non forested areas must, ideally, be removed.
#' @param res scalar. The resolution at which the estimation must be performed
#' @param f scalar. factor.
#' @param at numeric. densities to predict
#' @param normalize bool. Normalize on-the-fly if the dataset is not height normalized.
#'
#' @references Roussel, J.-R., Caspersen, J., Béland, M., Thomas, S., & Achim, A. (2017). Removing
#' bias from LiDAR-based estimates of canopy height: Accounting for the effects of pulse density and
#' footprint size. Remote Sensing of Environment, 198, 1–16. https://doi.org/10.1016/j.rse.2017.05.032
#'
#' @examples
#' library(lidR)
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' las = readLAS(LASfile)
#'
#' shapefile_dir <- system.file("extdata", "lake_polygons_UTM17.shp", package = "lidR")
#' lakes = sf::st_read(shapefile_dir)
#' las = classify_poi(las, LASWATER, roi = lakes)
#' las = filter_poi(las, Classification != LASWATER)
#'
#' map = hmax_error_map(las, 10)
#' q = quantile(map[], probs = 0.95, na.rm = TRUE)
#' map = terra::stretch(map, maxv = q)
#' plot(map)
#'
#' @export
hmax_error_map = function(las, res, at = c(1,2,5,10), f = 5, normalize = FALSE)
{
  if (methods::is(las, "LAS"))
  {
    las@data$res = res
    las@data$at = list(at)
    if (normalize) las = lidR::normalize_height(las, lidR::tin())
    ans = lidR::pixel_metrics(las, ~hmax_error_metric(X,Y,Z, res, at), res*f)
    ans = -ans
    return(ans)
  }
  else if (methods::is(las, "LAScatalog"))
  {
    options = list(raster_alignment = res*f)
    ans = lidR::catalog_map(las, hmax_error_map, res = res, f = f, at = at, .options = options)
    return(ans)
  }
  else
  {
    stop("Invalid input")
  }
}

#' Max height error estimation local
#'
#' This function evaluates, for a given pixel size the expected error made on the metrics hmax as a
#' function of the pulse density according to Roussel et al (2007) (see references). The metric hmax
#' is always underestimated because discrete laser pulses are statistically unlikely to hit the very
#' top of the canopy. The highest laser hit is necessarily below the theoretically highest point. How
#' much? It depends on the shape/complexity of the canopy and the sampling density. By evaluating the
#' shape of the canopy and using probability theory, it is possible to estimate the expected error
#' made as a function of the pulse density (see references).
#'
#' @return A list with 2 data.frame. The first contains the normalized histogram of probabilities
#' (figure 13) and the second contains the estimated errors (figure 14). See references.
#'
#' @param las LAS object from lidR
#' @param res scalar. The resolution at which the estimation must be performed
#'
#' @references Roussel, J.-R., Caspersen, J., Béland, M., Thomas, S., & Achim, A. (2017). Removing
#' bias from LiDAR-based estimates of canopy height: Accounting for the effects of pulse density and
#' footprint size. Remote Sensing of Environment, 198, 1–16. https://doi.org/10.1016/j.rse.2017.05.032
#'
#' @examples
#' library(lidR)
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' las = readLAS(LASfile)
#'
#' shapefile_dir <- system.file("extdata", "lake_polygons_UTM17.shp", package = "lidR")
#' lakes = sf::st_read(shapefile_dir)
#'
#' las = classify_poi(las, LASWATER, roi = lakes)
#' las = filter_poi(las, Classification != LASWATER)
#'
#' res = 20
#' ans = hmax_error_local(las, res)
#'
#' library(ggplot2)
#' # plot the histogram of probabilities (fig 13)
#' ggplot(ans$shape) +
#'  aes(x = height, y = probability) +
#'  geom_line() +
#'  coord_flip() +
#'  ylab("probability") +
#'  ylim(0,0.02) +
#'  xlab("Height") +
#'  theme_bw()
#'
#' # plot dependence of bias on pulse density (fig 14)
#' ggplot(ans$error) +
#'  aes(x = density, y = error) +
#'  geom_line() +
#'  scale_x_log10(breaks = c(1,2,3,4,5,10,20,30,50)) +
#'  scale_y_continuous(limits=c(-1.3,0), breaks =seq(-1.25,0,0.25)) +
#'  xlab("Pulse density") +
#'  ylab("Theorical error") +
#'  ggtitle(paste("Abacus computed for", res*res, "square meters")) +
#'  theme_bw()
#'
#' @export
hmax_error_local = function(las, res)
{
  X = Y = Z = . =  NULL
  if (methods::is(las, "LAS")) las = las@data

  bin  = 0.1
  A    = res*res
  zmax = ceiling((max(las$Z)+1)/bin)*bin   # zmax is higher than the highest Z value. See function canopy_histogram

  las$X = lidR:::f_grid(las$X, res, 0)
  las$Y = lidR:::f_grid(las$Y, res, 0)
  stat = las[, canopy_histogram(Z, zmax, bin), by = .(X,Y)]
  stat[, X := NULL]
  stat[, Y := NULL]

  # Compute the probabilities
  p = as.numeric(colSums(stat))
  h = as.numeric(names(stat))
  p = p/sum(p)

  # Compute the abacus for pulse densities from 1 to 50
  pulseDensities = seq(0.5, 50, 0.5)
  numOfPulses    = pulseDensities*A
  error          = sapply(numOfPulses, M, p, h)

  ans = list(
    shape = data.frame(height = h, probability = p),
    errors = data.frame(density = pulseDensities, error))

  return(ans)
}

#' Error metric
#'
#' Exported for internal use only
#'
#' @param X,Y,Z,res Internal use.
#' @param at numeric. Densities to estimate
#' @export
hmax_error_metric = function(X,Y,Z, res, at)
{
  res = res[1]
  at = at[1][[1]]
  data = data.table::data.table(X,Y,Z)
  ans = hmax_error_local(data, res)
  err = estimate_error(ans$errors, at)
  res = as.list(err$y)
  names(res) = err$x
  return(res)
}

estimate_error = function(abacus, density = c(0.5,1,2,3,4,6,8,10,15,20))
{
  stats::approx(abacus$density, abacus$error, density)
}
