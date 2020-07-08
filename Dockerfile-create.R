## Load / install required libraries
if (!(require(containerit))) {
  if (!(require(remotes))) install.packages("remotes") 
  remotes::install_github("o2r-project/containerit")
}
if (!(require(here))) install.packages("here") 

## Get the Dockerfile environment
my_env <- containerit::dockerfile(from = here::here("R/sceptic.R"))
## Print it nicely
cat(as.character(format(my_env)), sep = "\n")

## Write the Dockerfile to disk
write(my_env, file = here::here("Dockerfile"))