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
stan: $(standir)mod-pred.stan $(standir)mod.stan $(standir)mod-anthro.stan \
 $(standir)mod-me.stan
main: $(resdir)main/tcr.fst $(resdir)main/gmst2100.fst \
 $(resdir)main/gmst-pred.csv $(resdir)main/params.csv $(resdir)main/had-dev.csv
sensitivity: $(resdir)sensitivity/params-alt-gmst.csv $(resdir)sensitivity/tcr-alt-gmst.fst \
 $(resdir)sensitivity/params-me-gmst.csv $(resdir)sensitivity/tcr-me-gmst.fst \
 $(resdir)sensitivity/tcr-me-forcings.fst \
 $(resdir)sensitivity/tcr-eff1-forcings.fst $(resdir)sensitivity/tcr-eff2-forcings.fst \
 $(resdir)sensitivity/params-anthro.csv $(resdir)sensitivity/tcr-anthro.fst
 
clean:
	rm -f $(results_main) $(datdir)* $(rawdir)*

## Draw the Makefile DAG
## Requires: https://github.com/lindenb/makefile2graph
dag: makefile-dag.png
makefile-dag.png: Makefile
	make -Bnd all | make2graph | dot -Tpng -Gdpi=300 -o makefile-dag.png

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

#### b) Measurement error in GMST
results_me_gmst = $(resdir)sensitivity/params-me-gmst.csv $(resdir)sensitivity/tcr-me-gmst.fst
$(results_me_gmst) &: $(rdir)03-sensitivity-me-gmst.R $(standir)mod-me.stan $(datdir)climate.csv
	Rscript $<

#### c) Measurement error in forcings
results_me_gmst = $(resdir)sensitivity/tcr-me-forcings.fst
$(results_me_gmst) &: $(rdir)03-sensitivity-me-forcings.R $(standir)mod.stan $(datdir)climate.csv \
 $(datdir)df18.fst
	Rscript $<

#### d) Adjust forcing efficacies (Marvel et. al, 2016)
results_eff = $(resdir)sensitivity/tcr-eff1-forcings.fst $(resdir)sensitivity/tcr-eff2-forcings.fst
$(results_eff) &: $(rdir)03-sensitivity-eff.R $(standir)mod.stan $(datdir)climate.csv \
 $(rawdir)rcps.csv
	Rscript $<

#### e) Separate out anthropogenic forcings
results_anthro = $(resdir)sensitivity/params-anthro.csv $(resdir)sensitivity/tcr-anthro.fst
$(results_anthro) &: $(rdir)03-sensitivity-anthro.R $(standir)mod-anthro.stan $(datdir)climate.csv
	Rscript $<
	
	
## Helpers
.PHONY: all clean dag data stan main sensitivity
.DELETE_ON_ERROR:
.SECONDARY: