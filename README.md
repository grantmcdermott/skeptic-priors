# Sceptic priors and climate consensus

This repository contains code and data for my working paper, ["Sceptic priors and climate consensus"](https://grantmcdermott.com/papers/scepticpriors.pdf). 

> **Abstract:** How much evidence would it take to convince sceptics that they are wrong about climate change? I explore this question within a Bayesian framework. I consider a group of stylised sceptics and examine how these individuals update their beliefs in the face of current and continuing climate change. I find that available evidence in the form of instrumental climate data tends to overwhelm all but the most extreme priors. Most sceptics form updated beliefs about climate sensitivity that correspond closely to estimates from the scientific literature. However, belief convergence is a non-linear function of prior strength. It becomes increasingly difficult to convince the remaining pool of sceptics that they are wrong. I discuss the necessary conditions for consensus formation under Bayesian learning and show how apparent deviations from the Bayesian ideal still be accommodated within the same conceptual framework. I argue that a generalized Bayesian model thus provides a bridge between competing theories of climate scepticism as a social phenomenon.

Click on the green "Code" button above to clone or download the repo to your local computer. Alternately, click on the "fork" button at the very top right of the page to create an independent copy of the repo within your own GitHub account.

The scripts for running the analysis can be found in the `R/` sub-directory. Click on this sub-directory to see more details in the accompanying README file. However, first you need to make sure that you have completed the necessary software installation (below).
 
## Software requirements and reproducibility

Users have two options for recreating the environment needed to reproduce the analysis: i) Manual configuration, or ii) Docker.

### Manual configuration

#### Step 1. Install *R* (and RStudio)

All of the analysis is conducted in the *R* programming environment. *R* is free, open-source and available for download [**here**](https://www.r-project.org/). The code has been tested against R versions 3.5.0 - 3.6.3, although I would just as soon recommend using the latest version of R (4.0.2 at the time of writing).

*Optional:* I normally recommend running *R* in the RStudio IDE, which you can also download for free [**here**](https://www.rstudio.com/products/rstudio/download/). However, note that the most computationally-intensive model runs should be called directly from the terminal (e.g. using `Rscript`), since they employ a parallel forking process that, while very efficient, can cause problems if run through an IDE like RStudio. More details are provided in the relevant README files (e.g. [here](https://github.com/grantmcdermott/sceptic-priors/blob/master/R/Evidence/README.md)).

#### Step 2. Install JAGS

In addition, you will need to install JAGS ("Just Another Gibbs Sampler"), which is the underlying program used for running the Bayesian regressions. JAGS too is free and open-source, and is available for download [**here**](http://mcmc-jags.sourceforge.net/).

#### Step 3. Install *R* packages

Once *R* and JAGS are successfully set up on your system, you will need to install a number of external *R* packages. I have used [**renv**](https://rstudio.github.io/renv/) to create a lockfile that snapshops the necessary package versions. Installing the correct packages should thus be as simple as:

```r
#renv::init()    ## Only necessary if you didn't clone/open the repo as an RStudio project
renv::restore()  ## Enter "y" when prompted
```

*Optional:* The `extrafont` package is used to embed Fira Sans fonts in the figures. Please note that the Fira Sans font family must be installed separately on your system (e.g. from [here](https://fonts.google.com/specimen/Fira+Sans)) and also requires some minor setup before *R* recognizes it (instructions [here](https://github.com/wch/extrafont/blob/master/README.md)). However, you can also skip this setup if you want; the code is written in such a way that it will revert to *R*'s default Arial font if Fira Sans is unavailable.

### Docker

For those of you who don't feel like configuring a manual setup, I also provide a [Docker image](https://hub.docker.com/repository/docker/grantmcd/sceptic) that bundles all of the necessary elements. You'll also want to mount this repo on the running container so that you can access the actual data and scripts. For example, say that you've cloned this repo to `/home/yourname/sceptic-priors` on your computer. Then you would run

```sh
docker run -it --rm -v /home/yourname/sceptic-priors:/mnt/sceptic-priors grantmcd/sceptic /bin/bash
```

(This will take a few minutes to download the container image from DockerHub the first time you run it. But the container will be ready and waiting for immediate deployment on your system thereafter.)

Once the container is up and running, you should see something like:

```sh
root@53f880d5ea69:/payload# 
```

You can then navigate to the mounted volume and start running the scripts, e.g.

```
root@53f880d5ea69:/payload# cd /mnt/sceptic-priors/ 
root@53f880d5ea69:/mnt/sceptic-priors# Rscript --vanilla R/sceptic.R
```

Of course, you can also open up the interactive R shell (type `R`) instead of running the scripts. Exit and stop the Docker container by typing `exit`.

If you are unfamiliar with Docker, then I recommend starting [here](https://ropenscilabs.github.io/r-docker-tutorial/).

## Problems

If you have any trouble running the code, or find any errors, please file an issue on this repo and I'll look into it.

## License

The software code contained within this repository is made available under the [MIT license](http://opensource.org/licenses/mit-license.php). The data and figures are made available under the [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).
