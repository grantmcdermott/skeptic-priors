# Sceptic priors and climate policy

This repository contains code and data for replicating my working paper, ["Sceptic priors and climate policy"](https://drive.google.com/file/d/0B6AgOxtQA9dTcjRmZkNjMVhuVFU/view?usp=sharing). 

> **Abstract:** How much evidence would it take to convince climate sceptics that they are wrong about global warming? I explore this question within a Bayesian framework. I consider a group of stylised sceptics and examine how these individuals update their beliefs in the face of current and continuing climate change. I find that available evidence in the form of instrumental climate data already tends to overwhelm all but the most extreme priors. The resulting posterior distributions of climate sensitivity correspond closely to existing estimates from the scientific literature. The updated beliefs of most sceptics are thus consistent with a carbon price that is substantially greater than zero. However, belief convergence is a non-linear function of prior strength, so that it become increasingly difficult to convince the marginal sceptic. I conclude by discussing the general conditions for consensus formation under Bayesian learning, its relevance to our current policy impasse, and offer some remarks about finding common ground in the future.

Click on the "fork" button at the very top right of the page to create an independent copy of the repo within your own GitHub account. Alternately, click on the green "clone or download" button just below that to download the repo to your local computer.

The scripts for running the analysis can be found in the `R/` sub-directory. Click on this directory to see more details in the accompanying README file. However, first you need to make sure that you have completed the necessary software installation (below).

## Software requirements

### Step 1. Install *R* (and RStudio)

All of the analysis is conducted in the *R* programming environment. *R* is free, open-source and available for download [**here**](https://www.r-project.org/).  

*Optional:* I normally recommend running *R* in the RStudio IDE, which you can also download for free [**here**](https://www.rstudio.com/products/rstudio/download/). However, note that the most computationally-intensive model runs should be called directly from the terminal (e.g. using `Rscript`), since they employ a parallel forking process that, while very efficient, can cause problems if run through an IDE like RStudio. More details are provided in the relevant README files (e.g. [here](https://github.com/grantmcdermott/sceptic-priors/blob/master/R/Evidence/README.md)).

### Step 2. Install JAGS

In addition, you will need to install JAGS ("Just Another Gibbs Sampler"), which is the underlying program used for running the Bayesian regressions. JAGS too is free and open-source, and is available for download [**here**](http://mcmc-jags.sourceforge.net/).

### Step 3. Install *R* packages

Once *R* and JAGS are successfully set up on your system, you will need to install a number of external *R* packages. These are listed at the top of the `R/sceptic_funcs.R` script. However, a convenient way to ensure that you have the correct versions of all these packages is to simply run the following code chunk in your *R* console:

```r
if (!require("pacman")) install.packages("pacman")
pacman::p_install(c(LearnBayes, rjags, R2jags, dclone, snow, grid, gridExtra, tidyverse, devtools, hrbrthemes, ggridges, RColorBrewer, stargazer, xtable, pbapply, tictoc, extrafont, R.cache, here, RhpcBLASctl))
devtools::install_github("johnbaums/jagstools")
pacman::p_update()
```

The `extrafont` package is used to embed Fira Sans fonts in the figures. Please note that the Fira Sans font family must be installed separately on your system (e.g. from [here](https://fonts.google.com/specimen/Fira+Sans)) and also requires some minor setup before *R* recognizes it (instructions [here](https://github.com/wch/extrafont/blob/master/README.md)). However, you can also skip this setup if you want; the code is written in such a way that it will revert to *R*'s default Arial font if Fira Sans is unavailable.

## Problems

If you have any trouble running the code, or find any errors, please file an issue on this repo and I'll look into it.

## License

The software code contained within this repository is made available under the [MIT license](http://opensource.org/licenses/mit-license.php). The data and figures are made available under the [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).
