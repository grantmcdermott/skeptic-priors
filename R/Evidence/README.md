# Evidence needed for sceptic beliefs to converge with the mainstream

This directory contains code for estimating the amount of evidence that different sceptics need to converge with mainstream beliefs about climate change. As described in the paper, this is defined to have occurred once the mean posterior TCR (for a given prior) equals either 1.3 &deg;C or 1.5 &deg;C. Similar to the [recursive estimates](https://github.com/grantmcdermott/sceptic-priors/tree/master/R/Recursive) described elsewhere in the repository, this section of the analysis works by iterating over different priors, one year at a time. The main differences between the two sections are:

1. Whereas the recursive section focuses specifically on the four main sceptics defined in the paper ("Moderate Lukewarmer", "Strong Denier", etc.), the evidence section iterates over a more granular range.
2. Whereas the recursive section uses historical data, the data used as evidence in this section are all simulated. That is, simulated in the sense that they are based on parameters obtained from the noninformative Bayesian regression (including error term) rather than explicitly observed in the dataset. The reason for using simulated data should be clear from the paper: The more hardcore sceptics will only converge with the mainstream once additional data has been accumulated in the future.<sup>[1](#myfootnote1)</sup>

The entire code for executing this section of the analysis is contained with the R script, **`evidence.R`**. Please note that this script should *not* be run from within RStudio, but rather from an R instance running in the shell (e.g. via `$ Rscript evidence.R`). See [Warning about RStudio](#warning-about-rstudio) below for more details.

## Performance

This section is the most computationally-intensive part of the project and takes several hours to complete even if it is run on a server. Expect it to take considerably longer on local machines like regular laptops or desktop computers. (See [here](http://grantmcdermott.com/2017/05/30/rstudio-server-compute-engine/) for details on how to create your own server using Google Compute Engine, including a one-year free trial.) That being said, I have tried to optimise the code in several ways to mitigate the pain:

- The code is run in parallel and will exploit all available multicore resources. However, see [Warning about RStudio](#warning-about-rstudio) below.
- A progress bar and onscreen messages will give you a sense of how long you can expect to wait.
- Caching is automatically enabled via the [R.cache package](https://cran.r-project.org/web/packages/R.cache/index.html). So you can immediately resume from where you were if something unexpected happens (e.g. a problem like crashing or timeout), or after you have run the function once.
- Finally, if you don't wish to run everything yourself, then you can simply read in my previously saved results [here](https://github.com/grantmcdermott/sceptic-priors/blob/master/Results/Evidence/tcr-evidence.csv).

## Warning about RStudio

I strongly recommend that you run this section from the shell, rather than from within the RStudio IDE. Assuming you are already in the right directory, the the simplest way to do this is by opening up the shell and calling the `evidence.R` script as follows:

```
$ Rscript evidence.R
``` 

(where the `$` denotes a shell prompt that should be ignored.)

**Reason: ** The the code will automatically invoke *forking* for the parallel cluster if it is available (i.e. on Linux and Mac). This is considerably (3x) faster and less memory intensive than the alternative *parallel socket* approach. (See [here](https://raw.githack.com/uo-ec607/lectures/master/12-parallel/12-parallel.html#forking_vs_sockets) for background.) However, is well known that forking can cause problems when it is invoked from within a GUI or IDE like RStudio. And, indeed, the script will likely hang at random points in the parallel loop if you try to run it in RStudio. The shell approach circumvents this problem. An alternative is to manually change the `cl_type` object from "FORK" to "PSOCK" (i.e. the default for Windows) before you run the main `evid_func()` function (see +/- line 75). However, this will mean that the script will take even longer, possibly days, to complete.

--

<a name="myfootnote1"><sup>[1]</sup></a> Of course, it is possible to iterate over the observed data and then switch to simulated data only in cases where the historical record is not sufficient. However, the code becomes clunky and potentially confusing. (Do you iterate backwards from the present day for the historical data, but forwards for the simulated future data?) More importantly, the results are very similar regardless of which approach you take. Email me if you would like to see for yourself.
