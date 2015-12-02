# Sceptic priors and climate policy

This repository contains code and data for replicating my ["Sceptic priors and climate policy"](https://drive.google.com/file/d/0B6AgOxtQA9dTcjRmZkNjMVhuVFU/view?usp=sharing) paper.

The main file for running the analysis is "sceptic.R". This parent file calls several subsidiary files ("jags-loop.R", "noninf-loop.R", "scep_funcs.R", etc.), which will loop over the different sceptic prior types and climate scenarios, and also run the MCMC simulations for obtaining the posterior distributions of the Bayesian regressions.

## A brief note on performance (length of MCMC chains and parallelisation/clustering)
I have put all the MCMC simulations into a nested loop as that makes it easier for anyone to read, understand, and run the code in a concise manner. At the same time, this also means that the code takes a few minutes to run over all sceptic types (i.e. priors). You can obviously speed things up by reducing the length of the MCMC chains (`chain_length <- ...` on +/- line 32 in the "sceptic.R" parent file). A related point is that the MCMC simulations are run in parallel in JAGS to help speed things up. Depending on how many processors your computer has, you can also improve performance by increasing the number of clusters/chains on the next line (`n_chains <- ...`) .
