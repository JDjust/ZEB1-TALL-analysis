args0 <- commandArgs(FALSE)
filearg <- grep("^--file=", args0, value=TRUE)
root <- if (length(filearg)) dirname(normalizePath(sub("^--file=", "", filearg[1]), winslash="/")) else getwd()
while (!file.exists(file.path(root,"src","figure_engine","render_figure.py"))) {
  up <- dirname(root)
  if (identical(up,root)) stop("Cannot find repository root.")
  root <- up
}
source(file.path(root,"src","render.R"))
render_locked(root,"S1",panel="A")
