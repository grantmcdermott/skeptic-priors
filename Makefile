## Directory vars (usually only these need changing)
rawdir = data/raw/
datdir = data/
rdir = R/
standir = stan/
resdir = results/

## See note about new grouped targets method, i.e. replacing ":" with "&:"
## https://stackoverflow.com/a/59877127/4115816

## Headline build
all: data stan main sensitivity

data: $(datdir)climate.csv $(datdir)priors.csv $(datdir)df18.fst
stan: $(standir)mod-pred.stan $(standir)mod.stan $(standir)mod-anthro.stan
main: $(resdir)main/tcr.fst $(resdir)main/gmst2100.fst \
 $(resdir)main/gmst-pred.csv $(resdir)main/params.csv $(resdir)main/had-dev.csv
sensitivity: $(resdir)sensitivity/params-alt-gmst.csv $(resdir)sensitivity/tcr-alt-gmst.fst \
 $(resdir)sensitivity/params-anthro.csv $(resdir)sensitivity/tcr-anthro.fst
 
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

### Sensitivity analysis results

#### a) Alt GMST series
results_gmst_alt = $(resdir)sensitivity/params-alt-gmst.csv $(resdir)sensitivity/tcr-alt-gmst.fst
$(results_gmst_alt) &: $(rdir)03-sensitivity-alt-gmst.R $(standir)mod.stan $(datdir)climate.csv
	Rscript $<

#### e) Separate out anthropogenic forcings
results_anthro = $(resdir)sensitivity/params-anthro.csv $(resdir)sensitivity/tcr-anthro.fst
$(results_anthro) &: $(rdir)03-sensitivity-anthro.R $(standir)mod-anthro.stan $(datdir)climate.csv
	Rscript $<
	
	
## Helpers
.PHONY: all clean data
.DELETE_ON_ERROR:
.SECONDARY: