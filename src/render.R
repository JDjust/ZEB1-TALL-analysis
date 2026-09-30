# R controls the verified Python engine. No data acquisition is performed.
render_locked <- function(root, figure, panel=NULL, journal="canonical", out=NULL) {
  py <- Sys.getenv("ZEB1_PYTHON", unset="")
  if (!nzchar(py)) {
    candidates <- c("D:/Python/python.exe", unname(Sys.which("python")))
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    if (!length(candidates)) stop("Set ZEB1_PYTHON to an installed Python executable.")
    py <- candidates[1]
  }
  argv <- c(shQuote(file.path(root,"src","figure_engine","render_figure.py")), figure)
  if (!is.null(panel)) argv <- c(argv,"--panel",panel)
  if (journal!="canonical") argv <- c(argv,"--journal",shQuote(journal))
  if (!is.null(out)) argv <- c(argv,"--out",shQuote(out))
  status <- system2(py,argv)
  if (!identical(status,0L)) stop("Plot rendering failed, exit code ",status)
  invisible(status)
}
