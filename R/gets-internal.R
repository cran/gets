.onAttach <- function(libname, pkgname)
{

  mylib <- dirname(system.file(package = "gets"))

  txt <- c("\n",
    paste(sQuote("gets"), "version 0.9\n"),
    "\n",
    paste("An R package for general-to-specific (gets) modelling and"),
    paste("indicator saturation methods, see", sQuote("help(gets)"), "for details."),
    "\n",
    paste("CRAN website: https://cran.r-project.org/package=gets"),
    paste("For an introduction (PDF): http://www.sucarrat.net/R/gets"),
    "\n")

  if(interactive() || getOption("verbose")){
    packageStartupMessage(paste(strwrap(txt, indent = 2,
      exdent = 4), collapse = "\n"))
  }

}
