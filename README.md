# Sceptic priors and climate policy

This repository contains code and data for replicating my ["Sceptic priors and climate policy"](https://drive.google.com/file/d/0B6AgOxtQA9dTcjRmZkNjMVhuVFU/view?usp=sharing){:target="_blank"} paper. Click on the "fork" button at the very top right of the page to create an independent copy of the repo within your own GitHub account. Alternately, click on the green "clone or download" button just below that to download the repo to your local computer.

The main file for running the analysis is `sceptic.R`. This file will execute a nested loop, where the outer loop is over different priors types and the inner loop is over different climate scenarios. During each loop, the code will call several subsidiary scripts (e.g. `scep_funcs.R`, `jags-loop.R`, etc.) to run the Bayesian regressions, save the posterior results for later, and export them as figures or .tex files. Assuming that you have installed all of the necessary programs and packages (see below), you should thus be able to reproduce all of the primary results simply by running the `sceptic.R` parent file.

In addition, a number of supplementary regressions and simulations are described in the "Evidence", "Recursive", and "Robustness" folders. All of these supplementary exercises are similarly self-contained in the sense that they should execute fully upon running a single parent script. See the respective README files for details.

Lastly, the "Data" folder contains the main dataset and the code needed to construct it from scratch. Any other folder or files that I have not described in detail should hopefully be self-explanatory (e.g. TablesFigures).

## Requirements

All of the analysis is conducted in the *R* programming environment. *R* is free, open-source and available for download [here](https://www.r-project.org/){:target="_blank"}. You will also need download and install [JAGS](http://mcmc-jags.sourceforge.net/){:target="_blank"} in order to run the Bayesian regressions.

Once *R* and JAGS are set up on your system, please install the following *R* packages. All of the packages are available on CRAN except where noted.

```r
library(readr)
library(LearnBayes)
library(rjags)
library(dclone)
library(snow)
library(devtools)
library(jagstools) ## devtools::install_github("johnbaums/jagstools")
library(ggplot2)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(extrafont) ## See note below
library(stargazer)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(pbapply)
```

A brief note on fonts and figures: The figures in this paper are produced using the `ggplot2` package, but incorporating [Palatino Linotype](http://www.myfontfree.com/palatino-linotype-myfontfreecom126f31679.htm){:target="_blank"} (or [Open Sans](https://fonts.google.com/specimen/Open+Sans){:target="_blank"}) fonts to match the paper's overall style. Most users will likely have the Palatino TFFs installed on their systems already. However, you will still need to register them to your *R* instance using the `extrafont` package. See [here](https://github.com/wch/extrafont){:target="_blank"} for instructions. If that all sounds like too much work, don't worry: The figures will revert to the ggplot2 default, as long as the `extrafont` package has been installed and loaded.

## Performance
As mentioned above, much of the code is contained within a nested loop. The goal is to encourage reproducibility by making the code easy to read and understand, say nothing of the ability to execute everything in a concise manner. Having said that, it can take a few minutes to run over the full set of prior types and climate scenarios. On my system (quad core CPU with 16GB RAM), the main loop only takes around two minutes to run. (Users with older machines can speed things up by reducing the length of the MCMC chains: `chain_length <- ...` on +/- line 15 of the `sceptic.R` parent file.) However, some of secondary analyses contained in the "Recursive" and "Evidence" folders take considerably longer to run. See the respective README files in those folders for more details, but consider this fair warning. One final thing to note is that the MCMC simulations are run in parallel in JAGS to help speed things up. The code will automatically detect how many CPUs your system has available and adjust the number of parallel processes accordingly.
