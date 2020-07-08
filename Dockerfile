FROM rocker/r-ver:3.5.1
LABEL maintainer="sceptic"
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y git-core \
	jags \
	libcairo2-dev \
	libxml2-dev \
	make \
	pandoc \
	pandoc-citeproc
RUN ["install2.r", "abind", "assertthat", "backports", "bindr", "bindrcpp", "boot", "broom", "cellranger", "cli", "coda", "colorspace", "crayon", "dclone", "digest", "dplyr", "evaluate", "extrafont", "extrafontdb", "fansi", "forcats", "gdtools", "generics", "ggplot2", "ggridges", "glue", "gridExtra", "gtable", "haven", "here", "hms", "hrbrthemes", "htmltools", "httr", "jsonlite", "knitr", "labeling", "lattice", "lazyeval", "LearnBayes", "lubridate", "magrittr", "Matrix", "modelr", "munsell", "nlme", "pbapply", "pillar", "pkgconfig", "plyr", "purrr", "R.cache", "R.methodsS3", "R.oo", "R.utils", "R2jags", "R2WinBUGS", "R6", "RColorBrewer", "Rcpp", "readr", "readxl", "remotes", "RhpcBLASctl", "rjags", "rlang", "rmarkdown", "rprojroot", "rstudioapi", "Rttf2pt1", "rvest", "scales", "stargazer", "stringi", "stringr", "tibble", "tictoc", "tidyr", "tidyselect", "tidyverse", "utf8", "viridisLite", "withr", "xfun", "xml2", "xtable"]
RUN ["installGithub.r", "johnbaums/jagstools@ca5ee61e4a534b46b810b20aceb3ef7dc4448e2b"]
WORKDIR /payload/
CMD ["R"]
