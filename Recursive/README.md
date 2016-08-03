# Recursive regressions

This folder contains code for running the set of recursive regressions described in the paper. You can run the entire code for this section by opening the `sceptic-recursive.R` parent file, which will then call various subsidiary files as needed.

The idea of the recursive estimates is to start with a small sample of observations close to the present day and then iterate backwards one year at a time until your sample extends over the full common historical dataset (i.e. 1866-2005). During each iteration, we record the TCR that obtains from running our Bayesian regression on that particular sub-sample of the data. This recursive process is incorporated within a larger function that loops over all prior types. The results are saved at each stage of the iteration and finally exported to file (`./Recursive/Data/tcr-rec-historic.csv`) for later convenience.

Side note: The "Animation" folder contains a series of figures that have been exported at different stages of the main loop. These figures show the full posterior TCR densities for each prior type at that point in the recursive process. The figures are obviously not all shown in the paper. Rather, they are intended for presentations, where they can be combined into an animation to visualise updating behaviour in our Bayesian framework.

## Performance
The entire recursive loop takes about 25 minutes to complete on my system (quad core CPU with 16GB RAM). It may take considerably longer on older machines. If you don't wish to re-run the recursive regressions yourself, then you can simply read in the previously saved results using the command `read_csv("./Recursive/Data/tcr-rec-historic.csv")`; see also line +/-105 of the `sceptic-recursive.R` file.

## "Historic" vs "future" recursions
The paper focuses on historic recursive regressions, which begin close to the present day and then iterate backwards through time to the beginning of the historical dataset. The code is set up to run this version as the default. However, it is also possible to iterate forwards into the future (using simulated data) by changing a few lines in the code. I included this option mostly to produce animated figures that I have used in some presentations (see above). However, I would not recommend it for other users. Rather see the "Evidence" folder of the main repository for simulations and analysis of future climate change.
