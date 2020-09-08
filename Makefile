## Directory vars (usually only these need changing)
rawdir = data/raw/
datdir = data/
rdir = R/
standir = stan/
resdir = results/

## Headline build
all: $(datdir)climate.csv $(datdir)df18.fst $(datdir)priors.csv

clean:
	rm -f $(results_main) $(data) $(rawfiles)

## Helpers
.PHONY: all clean data
.DELETE_ON_ERROR:
.SECONDARY:

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

## Main run
#main = $(rdir)02-main.R
#$(main): $(data) $(standir)mod-pred.stan

## Results
#results_main = $(resdir)tcr.fst $(resdir)y-dev.fst $(resdir)y-pred.fst
#$(results_main): $(main)
#	Rscript $<
