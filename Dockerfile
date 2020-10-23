FROM rocker/r-ver:4.0.2

LABEL maintainer="Grant McDermott <grantmcd@uoregon.edu>"

## Add git, Make, graphviz, Pandoc, etc. and Python (and install scipy)
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y git-core \
  python3.8 \
  python3-pip \
	graphviz \
	pandoc \
	pandoc-citeproc \
	libxml2-dev \
	libicu-dev
RUN pip3 install scipy==1.5.1

WORKDIR sceptic-priors
COPY . .

# Install renv and R packages
ENV RENV_VERSION 0.12.0
RUN echo "options(renv.consent = TRUE)" >> .Rprofile
RUN echo "options(RETICULATE_MINICONDA_ENABLED = FALSE)" >> .Rprofile
RUN R -e "install.packages('remotes', repos = c(RSPM = 'https://packagemanager.rstudio.com/all/latest'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e "renv::restore(confirm = FALSE)"
RUN R -e "cmdstanr::install_cmdstan(cores = 2)"

# Extra: makefile2graph for Make DAG
RUN cd /tmp
RUN git clone https://github.com/lindenb/makefile2graph.git
RUN cd makefile2graph
RUN make
#RUN make install ## Not sure why but this is generating build errors
RUN cd /sceptic-priors

CMD ["bash"]
