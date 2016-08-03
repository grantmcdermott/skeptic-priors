# Evidence needed for sceptic beliefs to converge with the mainstream

This folder contains code for estimating the amount of evidence that different sceptics need to converge with mainstream beliefs about climate change. As described in the paper, this is defined to have occurred once the mean posterior TCR (for a given prior) equals either 1.3 &deg;C or 1.5 &deg;C. Similar to the recursive estimates described elsewhere in this repository, this section of the analysis works by iterating over different priors, one year at a time. The main differences between the two sections are:

1. Whereas the recursive section focuses specifically on the four main sceptics defined in the paper (Moderate Lukewarmer, Strong Denier, etc.), the evidence section iterates over a more granular range.
2. Whereas the recursive section uses historical data, the data used as evidence in this section are all simulated. That is, simulated in the sense that they are based on parameters obtained from the noninformative Bayesian regression (including error term) rather than explicitly observed in the dataset. The reason for using simulated data should be clear from the paper: The more hardcore sceptics will only converge with the mainstream once additional data has been accumulated in the future.<sup>[1](#myfootnote1)</sup>

The entire code for executing this section of the analysis is contained with the file, `evidence.R`.

## Performance
 On my system (quad core CPU with 16GB RAM), the main function takes nearly an hour to complete. It may take considerably longer on older machines. (The function does come with a progress bar, so you should at least get a sense of how long you can expect to wait.) If you don't wish to re-run the evidence estimates yourself, then you can simply read in the previously saved results using the command `read_csv("Evidence/Data/tcr-evidence.csv")`; see line +/-224 of the `evidence.R` file.

 --

<a name="myfootnote1"><sup>[1]</sup></a> Of course, it is possible to iterate over the observed data and then switch to simulated data only in cases where the historical record is not sufficient. However, the code becomes clunky and potentially confusing (due you iterate backwards over the historical data and then forwards for the simulated future data?). More importantly, the results are very similar regardless of which approach you take. Email me if you would like to see for yourself.
