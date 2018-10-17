# Recursive regressions

This directory contains code for running the set of recursive regressions described in the paper. You can run the entire code for this section by opening the `sceptic-recursive.R` parent file, which will then call various subsidiary files as needed.

The idea of the recursive estimates is to start with a small sample of observations close to the present day and then iterate backwards one year at a time until your sample extends over the full common historical dataset (i.e. 1866-2005). During each iteration, we record the TCR that obtains from running our Bayesian regression on that particular sub-sample of the data. This recursive process is incorporated within a larger function that loops over all prior types. The results are saved at each stage of the iteration and finally exported to file for later convenience.

## Performance
The entire recursive loop takes about 25 minutes to complete on my system (quad core CPU with 16GB RAM). It may take considerably longer on older machines. If you don't wish to re-run the recursive regressions yourself, then you can simply read in my previously saved results [here](https://github.com/grantmcdermott/sceptic-priors/blob/master/Results/Recursive/tcr-rec-historic.csv).

## "Historic" vs "future" recursions
The paper focuses on historic recursive regressions, which begin close to the present day and then iterate backwards through time to the beginning of the historical dataset. The code is set up to run this version as the default. However, it is also possible to iterate forwards into the future (using simulated data) by changing a few lines in the code. I included this option mostly to produce animated figures that I have used in some presentations. However, I would not recommend it for other users. Rather see the `R/Evidence` sibling directory for simulations and analysis of future climate change.
