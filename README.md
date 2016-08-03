# Sceptic priors and climate policy

This repository contains code and data for replicating my paper, ["Sceptic priors and climate policy"](https://drive.google.com/file/d/0B6AgOxtQA9dTcjRmZkNjMVhuVFU/view?usp=sharing).

The main file for running the analysis is `sceptic.R`. This parent file calls several subsidiary files (e.g. `jags-loop.R`, `noninf-loop.R`, `scep_funcs.R`, etc.) as it loops over the different sceptic prior types and climate scenarios described in the paper. During each loop, the code will do things like run a Bayesian regression (using MCMC) to obtain the posterior distributions of key parameters, and then save these or export them as figures and .tex files. Provided that you have installed all of the necessary programs and packages (see below), you should thus be able to replicate all of the main results by simply executing the aforementioned `sceptic.R` file.

In addition to these main results, a number of supplementary regressions and simulations are described in the "Evidence", "Recursive", and "Robustness" folders. All of these supplementary exercises are self-contained in the sense that they should also run completely upon executing a single parent script. See the respective README files for details.

Lastly, the "Data" folder contains the main dataset and details on how it was constructed. Any other folder or files that I have not described in detail here should hopefully be self-explanatory (e.g. TablesFigures).

## Requirements

You will need **JAGS** on your system in order to run the Bayesian regressions. See [here](http://mcmc-jags.sourceforge.net/) for download and installation instructions.

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

A note on fonts and figures: The figures in this paper are produced using `ggplot2`, but with the "Palatino Linotype" font family replacing the ggplot2 default font for stylistic consistency. This makes use of the `extrafont` package by Winston Chang. See the [package homepage](https://github.com/wch/extrafont) for installation instructions and then check to make sure that your computer actually has the Palatino font types available. Alternatively, simply change the preferred font type on line 30 of the "sceptic.R" parent file.
The figures in this repo are produced using the *R* package `ggplot2`, but with the "Palatino Linotype" font family replacing the ggplot2 default font. (This is to ensure stylistic consistency with the main PDF document.) Most users will probably have the Palatino Linotype TFFs installed on their systems already (you can download them  [here](http://www.myfontfree.com/palatino-linotype-myfontfreecom126f31679.htm) if not). However, you still need to register them to your *R* instance using the `extrafont` package. Full instructions on how to do that are available [here](https://github.com/wch/extrafont). If this all sounds like too much work, don't worry: The figures will simply employ the ggplot2 default, as long as the `extrafont` package has been installed and loaded.

## Performance
As described above, I have put most of the code into a nested loop. This will hopefully make the code easier for most people to read and understand, say nothing of the ability to execute it in a concise manner. At the same time, the code does take a few minutes to run over all the different prior types and climate scenarios. On my system (quad core CPU with 16GB RAM), the main loop takes around two minutes to run. (If you're running the code on an older system then you can speed things up by reducing the length of the MCMC chains: `chain_length <- ...` on +/- line 15 of the `sceptic.R` parent file.) However, some of secondary analyses contained in the "Recursive" and "Evidence" folders take considerably longer to run. See the respective README files in those folders for more details, but consider this fair warning. One final thing to note is that the MCMC simulations are run in parallel in JAGS to help speed things up. The code will automatically detect how many CPUs your system has available and adjust the number of parallel processes accordingly.
