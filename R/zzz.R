#' @noRd
.onLoad <- function(libname, pkgname) {
  # 尝试加载动态库（如果存在）
  # 在开发模式下可能还未编译，所以使用 try
  try(library.dynam("asmcjr", pkgname, libname), silent = TRUE)
  
  utils::globalVariables(c(
    "idealpt", "density", "Density", "stimulus",
    "x", "y", "group", "pch"
  ))
  
  op <- options()
  op.asmcjr <- list(
    asmcjr.plot.theme = "classic",
    asmcjr.verbose = FALSE
  )
  toset <- !(names(op.asmcjr) %in% names(op))
  if(any(toset)) options(op.asmcjr[toset])
  
  invisible()
}

#' @noRd
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage(
      "\n",
      "================================================\n",
      "  asmcjr v", utils::packageVersion("asmcjr"), "\n",
      "  Analyzing Spatial Models of Choice and Judgment\n",
      "  (2nd Edition, 2024)\n",
      "================================================\n",
      "\n",
      "Quick start:\n",
      "  - Browse vignettes: browseVignettes('asmcjr')\n",
      "  - Package help: help(package = 'asmcjr')\n",
      "  - Report issues: https://github.com/uniofessex/asmcjr/issues\n"
    )
    
    if (!requireNamespace("rjags", quietly = TRUE)) {
      packageStartupMessage(
        "\nNote: Bayesian functions require JAGS: https://mcmc-jags.sourceforge.io/\n"
      )
    }
  }
}

#' @noRd
.onUnload <- function(libpath) {
  try(library.dynam.unload("asmcjr", libpath), silent = TRUE)
}