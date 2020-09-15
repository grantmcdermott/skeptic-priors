## Directory vars (usually only these need changing)
rawdir = data/raw/
datdir = data/
rdir = R/
standir = stan/
resdir = results/

## See note about new grouped targets method, i.e. replacing ":" with "&:"
## https://stackoverflow.com/a/59877127/4115816

## Headline build
all: data stan main robustness

data: $(datdir)climate.csv $(datdir)priors.csv $(datdir)df18.fst
stan: $(standir)mod-pred.stan $(standir)mod.stan
main: $(resdir)main/tcr.fst $(resdir)main/gmst2100.fst \
 $(resdir)main/gmst-pred.csv $(resdir)main/params.csv $(resdir)main/had-dev.csv
robustness: $(resdir)robustness/params-alt-gmst.csv $(resdir)robustness/tcr-alt-gmst.csv
 
clean:
	rm -f $(results_main) $(datdir)* $(rawdir)*

## Raw Data
raw: $(rdir)00-data-raw.R
	Rscript $<
	rm Rplots.pdf

## Prep Data
$(datdir)climate.csv: $(rdir)01-data-prep.R $(rawdir)*
	Rscript $<
	rm Rplots.pdf

$(datdir)priors.csv: $(rdir)01-data-prep.R $(datdir)climate.csv
	Rscript $<
	rm Rplots.pdf

$(datdir)df18.fst: $(rdir)01-data-prep.R $(rawdir)df18.idlsave
	Rscript $<
	rm Rplots.pdf

## Results

### Main results
results_main = $(resdir)main/tcr.fst $(resdir)main/gmst2100.fst \
 $(resdir)main/gmst-pred.csv $(resdir)main/params.csv $(resdir)main/had-dev.csv
$(results_main) &: $(rdir)02-main.R $(standir)mod-pred.stan $(datdir)climate.csv
	Rscript $<

### Robustness checks/results_main

#### a) Alt GMST series
results_gmst_alt = $(resdir)robustness/params-alt-gmst.csv $(resdir)robustness/tcr-alt-gmst.csv
$(results_gmst_alt) &: $(rdir)03-robustness-alt-gmst.R $(standir)mod.stan $(datdir)climate.csv
	Rscript $<
	
## Helpers
.PHONY: all clean data
.DELETE_ON_ERROR:
.SECONDARY: