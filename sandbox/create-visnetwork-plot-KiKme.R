# Load the data ---
cvn <- readRDS("../KiKme-analysis/results/KiKme-grid.rds")
cvn_reduced <- CVN::strip_cvn(cvn)
class(cvn_reduced) <- "cvn"
# create the basic plots
visnetwork_plots <- CVN::visnetwork_cvn(cvn_reduced)



