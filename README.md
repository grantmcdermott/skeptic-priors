# Sceptic priors and climate policy

This repository contains code and data for replicating my ["Sceptic priors and climate policy"](https://drive.google.com/file/d/0B6AgOxtQA9dTcjRmZkNjMVhuVFU/view?usp=sharing) paper.

The main file for running the analysis is "sceptic.R". This parent file calls several subsidiary files ("jags-loop.R", "noninf-loop.R", "scep_funcs.R", etc.), which will loop over the different sceptic prior types and climate scenarios, and also run the MCMC simulations for obtaining the posterior distributions of the Bayesian regressions.

Please note that you will have to download and install *JAGS* on your computer in order to run the MCMC simulations. See [here](http://mcmc-jags.sourceforge.net/) for instructions.

A similar, albeit more minor point is that all of the figures are exported with Palatino fonts. This make use of the `extrafont` package by Winston Chang. See the [package homepage](https://github.com/wch/extrafont) for installation instructions and then check to make sure that your computer actually has the Palatino font types available. Alternatively, simply change the preferred font type on line 30 of the "sceptic.R" parent file.

## A brief note on performance (length of MCMC chains and parallelisation/clustering)
I have put all the MCMC simulations into a nested loop as that makes it easier for anyone to read and understand, as well as run the full code in a concise manner. At the same time, this means that the code takes a few minutes to run over all sceptic types (i.e. priors). You can speed things up by reducing the length of the MCMC chains: `chain_length <- ...` on +/- line 38 in the "sceptic.R" parent file. (Note that the MCMC simulations are already run in parallel in JAGS to help speed things up. The code should automatically adjust the number of parallel processes according how many cores your computer has available.)
