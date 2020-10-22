
<!-- README.md is generated from README.Rmd. Please edit that file -->

Sceptic priors and climate consensus
====================================

<!-- badges: start -->
<!-- badges: end -->

This repository contains code and data for my working paper, [“Sceptic
priors and climate
consensus”](https://grantmcdermott.com/papers/scepticpriors.pdf).

> **Abstract:** How much evidence would it take to convince sceptics
> that they are wrong about climate change? I explore this question
> within a Bayesian framework. I consider a group of stylised sceptics
> and examine how these individuals update their beliefs in the face of
> current and continuing climate change. I find that available evidence
> in the form of instrumental climate data tends to overwhelm all but
> the most extreme priors. Most sceptics form updated beliefs about
> climate sensitivity that correspond closely to estimates from the
> scientific literature. However, belief convergence is a non-linear
> function of prior strength. It becomes increasingly difficult to
> convince the remaining pool of sceptics that they are wrong. I discuss
> the necessary conditions for consensus formation under Bayesian
> learning and show how apparent deviations from the Bayesian ideal
> still be accommodated within the same conceptual framework. I argue
> that a generalized Bayesian model thus provides a bridge between
> competing theories of climate scepticism as a social phenomenon.

Click on the green “Code” button above to clone or download the repo to
your local computer. Alternately, click on the “fork” button at the very
top right of the page to create an independent copy of the repo within
your own GitHub account.

Reproducibility
---------------

I use [**Make**](https://www.gnu.org/software/make/) to automate the
entire project. Assuming that you have installed Make and have met all
of the other dependencies — see [below](#dependencies) — the **TL;DR**
version for reproducing everything is:

    ## Run these commands in the shell
    git clone git@github.com:grantmcdermott/sceptic-priors.git
    cd sceptic-priors
    make

You can also limit the build to subsections of the project by passing
Make a relevant meta target, e.g.

-   `make data` will construct the dataset
-   `make main` will run the main analysis
-   `make sensitivity` will run all of the sensitivity analyses etc.

See the [Makefile](Makefile) to get a sense of the options. The
associated DAG is [here](makefile-dag.png) (warning: it’s complicated).

Dependencies
------------

While the entire project can be reproduced via a single call to `make`,
users must first satisfy various software dependencies. There are two
options: i) Manual configuration, or ii) Docker.

### Manual configuration

#### Step 1. Install R and R libraries

With the exception of a tiny bit of Python code, all of the analysis is
conducted in the R programming environment. R is free, open-source and
available for download [**here**](https://www.r-project.org/). The code
has been tested against R version 4.0.2.

Once R is successfully set up on your system, you will need to install a
number of external R libraries. I have used
[**renv**](https://rstudio.github.io/renv/) to snapshot the project’s R
environment. To install all of the necessary R libraries, you need
simply run the following command from your R console:

    ## Run these commands in R
    # renv::init()    ## Only necessary if you didn't clone/open the repo as an RStudio project
    renv::restore()  ## Enter "y" when prompted

#### Step 2. Install CmdStan

The workhorse Bayesian regressions for this project are passed from R to
[**CmdStan**](https://mc-stan.org/users/interfaces/cmdstan). This
enables the Bayesian MCMC computation to complete much, much faster than
it would otherwise. While R and CmdStan are two separate programs, the
easiest way to install the latter is from the former. Assuming that you
have completed Step 1 above, run the following line from your R console:

    ## Run this command in R
    cmdstanr::install_cmdstan(cores = 2)

#### Step 3. Optional(ish)

-   As mentioned, I use a tiny bit of Python code to extract some [IDL
    ‘save’ data](https://pypi.org/project/IDLSave/) as part of the data
    prep process. Assuming that you have cloned my repo as-is and did
    not delete any of the data files, Make will automatically skip this
    section and the Python requirement will be moot. Failing that,
    however, you will need **SciPy**’s [File
    IO](https://docs.scipy.org/doc/scipy/reference/tutorial/io.html)
    module. Regular Python users will almost certainly have SciPy
    installed on their system already. If not, you can install it
    yourself, e.g. with PyPi or Conda. I strongly recommend that R users
    go through **reticulate** (see
    [here](https://rstudio.github.io/reticulate/articles/python_packages.html)).
    FWIW, I am using Python 3.8.5 and SciPy 1.5.1 at the time of
    writing.

-   The `extrafont` package is used to embed [Fira
    Sans](https://fonts.google.com/specimen/Fira+Sans) fonts in the
    figures. Please note that the Fira Sans font family must be
    installed separately on your system and also requires some minor
    setup before R recognizes it (instructions
    [here](https://github.com/wch/extrafont/blob/master/README.md)).
    However, you can also skip this setup if you want; R’s default Arial
    fonts will be used if Fira Sans is unavailable.

### Docker

For those of you who don’t feel like configuring a manual setup, I also
provide a [Docker
image](https://hub.docker.com/repository/docker/grantmcd/sceptic) that
bundles all of the necessary elements. You’ll also want to mount this
repo on the running container so that you can access the actual scripts.
For example, say that you’ve cloned this repo to
`/home/yourname/sceptic-priors` on your computer. Then you would run

    docker run -it --rm -v /home/yourname/sceptic-priors:/mnt/sceptic-priors grantmcd/sceptic /bin/bash

(This will take a few minutes to download the container image from
DockerHub the first time you run it. But the container will be ready and
waiting for immediate deployment on your system thereafter.)

Once the container is up and running, you should see something like:

    root@53f880d5ea69:/payload# 

You can then navigate to the mounted volume and start running the
scripts, e.g.

    root@53f880d5ea69:/payload# cd /mnt/sceptic-priors/ 
    root@53f880d5ea69:/mnt/sceptic-priors# Rscript --vanilla R/sceptic.R

Of course, you can also open up the interactive R shell (type `R`)
instead of running the scripts. Exit and stop the Docker container by
typing `exit`.

If you are unfamiliar with Docker, then I recommend starting
[here](https://ropenscilabs.github.io/r-docker-tutorial/).

Performance
-----------

The code has been refactored to run all the Bayesian computation through
**CmdStan**. This has yielded *considerable* speed gains, to the point
that the entire analysis can be completed in under 30 minutes on my
[Dell Precision
5530](https://wiki.archlinux.org/index.php?title=Dell_Precision_5530)
laptop.<sup id="a1">[1](#f1)</sup> The table below provides a detailed
performance record for the different model runs. Note that I am
excluding the data preparation and paper production steps, but each of
these only takes a few seconds.

| Run                                         | File                             | Time (sec) | Cores used | RAM | OS         | Architecture                  |
|:--------------------------------------------|:---------------------------------|-----------:|-----------:|----:|:-----------|:------------------------------|
| Main                                        | `R/02-main.R`                    |      29.69 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |
| Recursive                                   | `R/03-recursive.R`               |     163.82 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |
| Evidence                                    | `R/04-evidence.R`                |     497.80 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |
| Alternative GMST series (sensitivity)       | `R/05-sensitivity-alt-gmst.R`    |       5.65 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |
| Anthrogenic forcings separate (sensitivity) | `R/05-sensitivity-anthro.R`      |       8.92 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |
| Adjusted forcings efficacies (sensitivity)  | `R/05-sensitivity-eff.R`         |     461.60 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |
| Measurement error in forcings (sensitivity) | `R/05-sensitivity-me-forcings.R` |     442.32 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |
| Measurement error in GMST (sensitivity)     | `R/05-sensitivity-me-gmst.R`     |       7.41 |         12 |  31 | Arch Linux | x86\_64-pc-linux-gnu (64-bit) |

Problems
--------

If you have any trouble running the code, or find any errors, please
file an issue on this repo and I’ll look into it.

License
-------

The software code contained within this repository is made available
under the [MIT license](http://opensource.org/licenses/mit-license.php).
The data and figures are made available under the [Creative Commons
Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).

<sup><b id="f1">1</b></sup> Previous versions of the code took a whole
day to complete on a cloud server. So you’ll understand my being pleased
by this improvement. [↩](#a1)
