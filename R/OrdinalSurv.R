#' Create a OrdinalSurv Object
#' @param id Observation subject's ID.
#' @param time Observation time.
#' @param ord Ordinal observation of event count
#' @param lower (Optional) lower cut points for ordinal observations
#' @param upper (Optional) upper cut points for ordinal observations
#' @return An object of S3 class \code{"OrdinalSurv"}.
#' @export
OrdinalSurv = function(id, time, ord, lower = NULL, upper = NULL) {
  os = list(osDF = data.frame(id = id, time = time, ord = ord))
  if (!is.null(lower) & !is.null(upper)) {
    os$osDF$lower = lower
    os$osDF$upper = upper
  }
  class(os) = "OrdinalSurv"
  os
}
