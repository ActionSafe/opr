#
# This is a Plumber API. You can run the API by clicking
# the 'Run API' button above.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#

library(plumber)
library(ggplot2)
library(dplyr)

#* @apiTitle Plumber Example API
#* @apiDescription Plumber example description.

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
function(msg = "") {
    list(msg = paste0("The message is: '", msg, "'"))
}

#* Plot a histogram
#* @serializer png list(width = 800, height = 600)
#* @post /plot
function(xvar = "Bill Length (mm)", yvar = "Bill Depth (mm)",
         species = c("Adelie", "Gentoo", "Chinstrap"),
         by_species = TRUE, show_margins = FALSE, smooth = TRUE) {
  subsetted = df |> filter(Species %in% species)
  p <- ggplot(subsetted, aes(x = !!sym(xvar), y = !!sym(yvar))) + list(
    theme(legend.position = "bottom"),
    if (by_species) aes(color = Species),
    geom_point(),
    theme_bw(base_size = 16),
    if (smooth) geom_smooth()
  )
  if (show_margins) {
    margin_type <- if (by_species) "density" else "histogram"
    p <- ggExtra::ggMarginal(p, type = margin_type, margins = "both",
                             size = 8, groupColour = by_species, groupFill = by_species)
  }
  print(p)
}

#* Return the sum of two numbers
#* @param a The first number to add
#* @param b The second number to add
#* @post /sum
function(a, b) {
    as.numeric(a) + as.numeric(b)
}

# Programmatically alter your API
#* @plumber
function(pr) {
    pr %>%
        # Overwrite the default serializer to return unboxed JSON
        pr_set_serializer(serializer_unboxed_json())
}
