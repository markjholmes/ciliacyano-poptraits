# Ciliates and cyanobacteria: populations and traits

Contains code and data necessary to generate figures associated with manuscript "Environment- and system-specific interactions between population and trait dynamics".

- 00-data.r : code to load in the raw data and transform into more useful forms including performing the PCA to aggregate the traits
- 01-fit-models.r : fitting the different models described in text i.e. growth ~ f(N, T), trait change ~ g(N, T)
- 02-bootstrapping.r : running the PBLR test
- 03-model-info.r : combining and summarising the test results
- 04-plots.r : create the figures
- data/ : data storage
- figures/ : (initially doesn't exist) contains figure outputs