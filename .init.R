# .init.R
# Functions to initialize this Exercise session
# Boris Steipe
# ====================================================================
file.copy(list.files("./assets", full.names=TRUE), tempdir())

panelViewer <- getOption("viewer")

ABC.fig <- function(name) {
  panelViewer(file.path(tempdir(), name), height=580)
}

file.edit("R_Exercise-Graphs.R")

# [End]
