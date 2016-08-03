# Robustness checks

This folder contains a secondary set of regressions and simulations, aimed at establishing the robustness of the main results. In particular:

1. The "cw2014-gistemp" sub-folder contains code for running the noninformative Bayesian regression described in the paper, except now using two alternate sets of GMST data (CW2014 and GISTEMP) instead of the dafault HadCRUT4 series. As described in the paper, this is primarily motivated by the fact that the HadCRUT4 series suffers from slight biases in global coverage.

2. The "MeasError" sub-folder contains code for incorporating knowledge about the measurement error (ME) of GMST data into the Bayesian regressions. This approach makes use of the fact that the main GMST products (HadCRUT4, CW2014 and GISTEMP) all come with historical estimates of ME, as well as the fact that specifying the ME is easily done within the nested structure of the Bayesian regression framework.
