# Robustness checks

This directory contains a set of alternative regression specifications and simulations, aimed at establishing the robustness of the main results. In particular:

1. The `CW2014-GISTEMP/` sub-directory contains code for running a noninformative Bayesian regression, except now using two alternate sets of GMST data (CW2014 and GISTEMP) instead of the default HadCRUT4 series. As described in the paper, this is primarily motivated by the fact that the HadCRUT4 series suffers from slight biases in global coverage.

2. The `MeasError/` sub-directory contains code for incorporating knowledge about the measurement error (ME) of GMST data into the Bayesian regressions. This approach makes use of the fact that the main GMST products (HadCRUT4, CW2014 and GISTEMP) all come with historical estimates of ME, as well as the fact that specifying the ME is easily done within the nested structure of the Bayesian regression framework.

3. The `Marvel/` sub-directory contains code for adjusting the efficacies of certain radiative forcings as per [Marvel *et al.* (2016)](http://dx.doi.org/10.1038/nclimate2888). Since Marvel *et al.* provide distributional information about their efficacy parameter estimates, two approaches are considered. The first approach takes the mean efficacy estimates as given and ignores all uncertainty. The second approach tries to account for uncertainty by sampling from the full efficacy parameter distributions.

You can call all three sub-routines directly from the **`robustness.R`** parent file.
