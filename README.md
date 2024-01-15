# Nascent Chromatin Capture Proteomics

Maintainer: Marek Gierlinski (M.Gierlinski@dundee.ac.uk)

Software to accompany manuscript Bandau et al. (2024) "RNA polymerase II has distinct roles in configuring nascent and steady state chromatin".

## Usage

We suggest using RStudio. Start in the top project directory. The first step is to create environment using 'renv':

```
install.packages("renv")
renv::restore()
```

This will install all necessary packages. Run the `targets` pipeline.

```
targets::tar_make()
```

This will analyse data, create figures (in directory `fig`) and a report (in directory `doc`). 
