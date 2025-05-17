# set cran mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))
install_if_missing <- function(package) {
  packages <- installed.packages()[, 1]
  if (!is.element("package", packages)) {
    install.packages(package, dependencies = TRUE)
  }
}
install_if_missing("rmarkdown")
install_if_missing("R.matlab")
install_if_missing("viridis")
install_if_missing("ggplot2")
install_if_missing("latex2exp")
install_if_missing("MASS")
install_if_missing("Matrix")
install_if_missing("openxlsx")
install_if_missing("mclust")
install_if_missing("clue")
install_if_missing("aricode")
install_if_missing("cluster")
install_if_missing("proxy")
install_if_missing("openxlsx")
install_if_missing("rio")
install_if_missing("tidyr")
install_if_missing("dplyr")
install_if_missing("kableExtra")
install_if_missing("ggpubr")
install_if_missing("gridExtra")
install_if_missing("magick")
install_if_missing("GFA")
