



# internal functions ----------------------------------------------------------



#' Check if number(s) is/are integers. In contrast to is.integer, is_wholenumber
#' does not check the class but accepts all numbers that are integers with reasonable
#' precision.
#' @param x Number to test
#' @param tol tolerance for testing
#' @NoRd
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol




# useful functions --------------------------------------------------------




#' @title Feature Plot
#' @description Highligh gene expression data in a 2D-embedding (UMAP, tSNE< etc.).
#' @param embedding A matrix/data.frame/tibble (or anything that can be reasonably
#' converted with \code{data.frame(embedding)}) with exactly two columns.
#' If colnames are missing, the axis will be named "Dim1" and "Dim2".
#' @param expression Numeric vector with expression values of the gene of interest.
#' one element per cell. Usually normalized, recommended normalization is
#' \code{k/s/mean(1/s)}, where \code{k} are UMI counts for the gene of interest
#' and \code{s} are totalUMI of the cell (aka library size).
#' @param legend_name Text displayed above the legend. Most commonly the name
#' of the displayed gene.
#' @return A \code{ggplot2} object storing a colored scatter plot.
#' @details This function allows little customization on purpose, because it
#' bundles geoms, themes and settings that I found important for
#' visualizing gene expression in scRNAseq data:
#'
#' \itemize{
#'  \item coord_fixed, to avoid distortion of embeddings
#'  \item geom_point with size=.4, to ameliorate overplotting
#'  \item No axis ticks or background grid, because distances and axis units
#'  in embeddings do not
#'  carry meaning for most dimensionality reduction techniques.
#'  \item Intensity-coded color scales (viridis) displayed with
#'  log2-transformation. Makes visualization independent of colorblindness
#'  and appropriate for gene expression data (which is usually Log Normal
#'  distributed).
#'  \item Color scale breaks are displayed as 'closed interval', i.e.
#'  \code{max(expression)} and \code{min(expression)} are the most extreme
#'  breaks. Rounding makes them human-readable. This functionality is provided
#'  by \link[scUtils]{closed_breaks_log} and \link[scUtils]{closed_labels}.
#'       }
#'
#' If you insist on customizing, think of this function as a great starting point, you can simply
#' copy-paste the code after typing \code{feat} into your
#' console.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # expression goes from 0 to 22:
#'  set.seed(100)
#'  feat(matrix(rnorm(2000, c(.1, 3)), ncol=2), rpois(1000, c(.1, 11)))
#'  # expression goes from 2 to 52:
#'  set.seed(100)
#'  feat(matrix(rnorm(2000, c(.1, 3)), ncol=2), rpois(1000, c(10, 31)))
#'  }
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}}
#'  \code{\link[scUtils]{closed_labels}}
#'  \code{\link[scUtils]{closed_breaks_log}}
#' @rdname feat
#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols
#' @importFrom ggplot2 aes
feat <- function(embedding, expression, legend_name="Expression") {
  assertthat::assert_that(is.vector(expression), is.numeric(expression))
  assertthat::assert_that(is.numeric(embedding[,1, drop=TRUE]),
                          is.numeric(embedding[,2, drop=TRUE]),
                          ncol(embedding)==2)
  assertthat::assert_that(nrow(embedding) == length(expression))
  assertthat::assert_that(assertthat::is.string(legend_name))
  assertthat::assert_that(all(expression >= 0))
  # save axis_names for later:
  axis_names <- if (!is.null(colnames(embedding))){colnames(embedding)}else{
    c("Dim1", "Dim2")}
  colnames(embedding) <- c("u1", "u2")
  u1 <- u2 <- NULL # avoids "no visible binding" warning when building package
  # convert to data.frame (e.g. to handle matrices, etc.)
  embedding <- data.frame(embedding)

  # avoid zeros as it causes errors with log2 breaks. Replacing with a 10th of
  # expression's non-zero minimum seems reasonable:
  has_zeros <- any(expression == 0)
  expression[expression==0] <- min(expression[expression!=0])/10

  ggplot2::ggplot(data = dplyr::bind_cols(embedding, expression=expression),
                  mapping = ggplot2::aes(u1, u2, col=expression)) +
      ggplot2::geom_point(size=.4) + ggplot2::coord_fixed() +
    ggplot2::xlab(axis_names[1]) + ggplot2::ylab(axis_names[2])+
    viridis::scale_color_viridis(
      trans="log2",
      breaks = function(lims) closed_breaks_log(lims, base=2),
      labels=function(br) closed_labels(br, min_is_zero = has_zeros),
      na.value=viridisLite::viridis(1),
      name = legend_name)

}



#' @title Closed Breaks for log scales
#' @description Same algorithm as in \link[scales]{breaks_log}, but
#' forces inclusion of upper and lower limits (displaying the closed interval).
#' Including boundaries is particularly useful for ggplot2's color/fill, as it
#' emphasizes the meaning of maximal/minimal color intensities, see examples.
#' @param lims Vector with lower and upper limits of the data that you want
#' breaks for.
#' @param n Number of breaks to aim for. The actual number will be what the
#' algorithm can accomodate, see \code{?scales::breaks_log()}. Default: 5
#' @param base The base of logarithm to use. Default: 10
#' @return Numeric vector with breaks.
#' @details Code adapted from scales::breaks_log, original code is
#' [licensed by Hadley Wickham](https://cran.r-project.org/web/packages/scales/LICENSE),
#' changes are licensed by Felix Frauhamer.
#' Original license is MIT:
#'
#' \itemize{
#'   \item YEAR: 2010-2016
#'   \item COPYRIGHT HOLDER: Hadley Wickham}
#'
#'
#' The scUtils pacakges uses this function to color by gene expression,
#' where the maximal expression gives valuable
#' intuition for a gene's overall expression strength. For x- or y-axis, I still
#' recommend breaks_log from the scales package.
#' @examples
#' \dontrun{
#' if(interactive()){
#' # closed breaks include maximum, breaks_log do not:
#' closed_breaks_log(lims = c(.01, 977.1))
#' scales::breaks_log()(c(.01, 977.1))
#' # closed_breaks_log intends to replace outer breaks with actual limits:
#' round(scales::breaks_log(base=2)(c(.01, 977.1)), 3)
#' closed_breaks_log(lims = c(.01, 977.1), base=2)
#'  }
#' }
#' @seealso
#'  \code{\link[assertthat]{assert_that}}
#' @rdname closed_breaks_log
#' @export
#' @importFrom assertthat assert_that
closed_breaks_log <- function (lims, n=5, base=10)
{
  raw_rng <- suppressWarnings(range(lims, na.rm = TRUE))
  if (any(!is.finite(raw_rng))) {
    return(numeric())
  }
  assertthat::assert_that(min(lims) > 0) # added this line
  rng <- log(raw_rng, base = base)
  min <- floor(rng[1])
  max <- ceiling(rng[2])
  if (max == min)
    return(base^min)
  by <- floor((max - min)/n) + 1
  breaks <- base^seq(min, max, by = by)
  # remove breaks too close to original boundaries:
  breaks <- breaks[breaks > base * base^rng[1] &  breaks < base^rng[2]/base ]
  # specifically add original min and max:
  breaks <- c(base^rng[1], breaks, base^rng[2])
  relevant_breaks <- base^rng[1] <= breaks & breaks <= base^rng[2]
  if (sum(relevant_breaks) >= (n - 2)) {
    return(breaks)}
  while (by > 1) {
    by <- by - 1
    breaks <- base^seq(min, max, by = by)
    relevant_breaks <- base^rng[1] <= breaks & breaks <=
      base^rng[2]
    if (sum(relevant_breaks) >= (n - 2))
      return(breaks)
  }

  return(log_sub_breaks(rng, n = n, base = base))
}



#' @author Thierry Onkelinx, \email{thierry.onkelinx@inbo.be}
#' @noRd
#' @description Code adopted from the R package "scales", version 1.1.1.9000.
#' License of the scales package is MIT with:
#' YEAR: 2010-2016
#' COPYRIGHT HOLDER: Hadley Wickham
#' @param rng The range.
#' @param n Number of breaks, Default: 5
#' @param base base of logarithm, Default: 10
#' @return Vector of breaks
#' @details Internal function.
#' @examples
#' @seealso
#'  \code{\link[scales]{breaks_extended}}
#' @importFrom scales extended_breaks
log_sub_breaks <- function (rng, n = 5, base = 10)
{
  min <- floor(rng[1])
  max <- ceiling(rng[2])
  if (base <= 2) {
    return(base^(min:max))
  }
  steps <- 1
  delta <- function(x) {
    min(diff(log(sort(c(x, steps, base)), base = base)))
  }
  candidate <- seq_len(base)
  candidate <- candidate[1 < candidate & candidate < base]
  while (length(candidate)) {
    best <- which.max(vapply(candidate, delta, 0))
    steps <- c(steps, candidate[best])
    candidate <- candidate[-best]
    breaks <- as.vector(outer(base^seq(min, max), steps))
    relevant_breaks <- base^rng[1] <= breaks & breaks <=
      base^rng[2]
    if (sum(relevant_breaks) >= (n - 2)) {
      break
    }
  }
  if (sum(relevant_breaks) >= (n - 2)) {
    breaks <- sort(breaks)
    lower_end <- pmax(min(which(base^rng[1] <= breaks)) -
                        1, 1)
    upper_end <- pmin(max(which(breaks <= base^rng[2])) +
                        1, length(breaks))
    breaks[lower_end:upper_end]
  }
  else {
    scales::extended_breaks(n = n)(base^rng)
  }
}





#' @title Human-readable labels for closed breaks
#' @description Complements the closed_breaks_log function.
#' @param x Vector of breaks for which to produce labels.
#' Typically, this is the output of \code{closed_breaks_log}.
#' @param min_is_zero Should the smallest break be
#' displayed as zero (TRUE) or as the actual value (FALSE). Default: FALSE
#' @return Character vector with labels, used by \code{feat} function.
#' @details This is a helper for the \code{feat} function.
#' \code{feat} replaces numeric zeros with the next-smallest expression value
#' to avoid taking the logarithm of zero. \code{min_is_zero} can be used to
#' display the lowest break of the color scale as zero in these cases.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # human readable output:
#'  closed_labels(c(.001111,.122, 0.5, 10, 100, 1800))
#'  }
#' }
#' @seealso
#'  \code{\link[scales]{label_scientific}}
#'  \code{\link[scales]{label_number_auto}}
#' @rdname closed_labels
#' @export
#' @importFrom dplyr case_when
#' @importFrom scales scientific
closed_labels <- function(x, min_is_zero = FALSE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr must be installed for this functionality.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("scales must be installed for this functionality.")
  }
  x <- dplyr::case_when(
    x == min(x) & min_is_zero ~ as.character(0), # for feat function
    x == 0 ~ as.character(0),
    # I want maximum to be labelled more precicely than others:
    x == max(x) & base::abs(x) >= 1000 ~ scales::scientific(x,  digits = 2),
    base::abs(x) < .01 | base::abs(x) >= 1000 ~ scales::scientific(x,  digits = 0),
    x > 2 ~ as.character(round(x, digits=1)),  #
    TRUE ~ as.character(round(x, 2)))
  # remove leading zero because it looks nicer:
  return(base::gsub("^0.", ".", x))
}







# quick tests -------------------------------------------------------------
# see if breaks and labels look pleasant
# check_pleasentness <- function(lim1, lim2, base=10){
#   ggplot2::ggplot(data= data.frame(x=1:15, y =seq(lim1, lim2, length.out = 15)),
#          ggplot2::aes(x, x, col=y))+ggplot2::geom_point(size=5) +
#   viridis::scale_color_viridis(trans="log2",
#                        breaks=function(lims) {print(lims)
#                          closed_breaks_log(lims, n=5, base)},
#                        labels = closed_labels)
#
# }
# # compare to scales::label_number_auto
# check_pleasentness_auto <- function(lim1, lim2, base=10){
#   ggplot2::ggplot(data= data.frame(x=1:15, y =seq(lim1, lim2, length.out = 15)),
#                   ggplot2::aes(x, x, col=y))+ggplot2::geom_point(size=5) +
#     viridis::scale_color_viridis(trans="log2",
#                                  breaks=function(lims) {print(lims)
#                                    closed_breaks_log(lims, n=5, base)},
#                                  labels = scales::label_number_auto())
#
# }
# scUtils:::check_pleasentness(.244, 977.223,2)
# scUtils:::check_pleasentness(.244, 101.223,2)
# scUtils:::check_pleasentness(.9, 200.1, 2)
# scUtils:::check_pleasentness(.12444, 977.223, base=2)
# scUtils:::check_pleasentness(.1,.244, 2)







