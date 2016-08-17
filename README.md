# Sceptic priors and climate policy

This repository contains code and data for replicating my paper, ["Sceptic priors and climate policy"](https://drive.google.com/file/d/0B6AgOxtQA9dTcjRmZkNjMVhuVFU/view?usp=sharing).

The main file for running the analysis is `sceptic.R`. This parent file calls several subsidiary files (e.g. `scep_funcs.R`, `jags-loop.R`, etc.) as it loops over the different prior types and climate scenarios described in the paper. During each loop, the code will do things like run a Bayesian regression (using MCMC) to obtain the posterior distributions of key parameters, save these results for later, or export them as figures and .tex files. Assuming that you have installed all of the necessary programs and packages (see below), you should be thus able to replicate all of the main results by simply executing the `sceptic.R` parent file.

In addition to the main results, a number of supplementary regressions and simulations are described in the "Evidence", "Recursive", and "Robustness" folders. All of these supplementary exercises are similarly self-contained in the sense that they should also run completely upon executing a single parent script. See the respective README files for details.

Lastly, the "Data" folder contains the main dataset and details on how it was constructed. Any other folder or files that I have not described in detail here should hopefully be self-explanatory (e.g. TablesFigures).

## Requirements

All of the analysis is conducted in the *R* programming environment. *R* is free, open-source and available for download [here](https://www.r-project.org/). You will also need download and install [JAGS](http://mcmc-jags.sourceforge.net/) in order to run the Bayesian regressions. Click the link for instructions.

In addition to JAGS, please make sure the following *R* packages are installed on your system. All of the packages are available on CRAN except where noted.

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

A note on fonts and figures: The figures in this paper are produced using the `ggplot2` package, but incorporating a [Palatino Linotype](http://www.myfontfree.com/palatino-linotype-myfontfreecom126f31679.htm) font to match the paper's overall look. Most users will likely have the Palatino TFFs installed on their systems already. However, you will need to register them to your *R* instance using the `extrafont` package. See [here](https://github.com/wch/extrafont) for instructions. If that all sounds like too much work, don't worry: The figures will simply revert to the ggplot2 default, provided the `extrafont` package has been installed and loaded.

## Performance
As mentioned above, I have written most of the code in a nested loop. The goal is to make the code easier to read and understand, say nothing of the ability to execute in a concise manner. Having said that, it can take a few minutes to run over the full set of prior types and climate scenarios. On my system (quad core CPU with 16GB RAM), the main loop only takes around two minutes to run. (Users with older machines can speed things up by reducing the length of the MCMC chains: `chain_length <- ...` on +/- line 15 of the `sceptic.R` parent file.) However, some of secondary analyses contained in the "Recursive" and "Evidence" folders take considerably longer to run. See the respective README files in those folders for more details, but consider this fair warning. One final thing to note is that the MCMC simulations are run in parallel in JAGS to help speed things up. The code will automatically detect how many CPUs your system has available and adjust the number of parallel processes accordingly.
