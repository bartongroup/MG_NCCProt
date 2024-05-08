# Nascent Chromatin Capture Proteomics

Maintainer: Marek Gierlinski (M.Gierlinski@dundee.ac.uk)

Collaborators: Susanne Bandau, Constance Alabert

Software to accompany manuscript Bandau et al. (2024) "RNA polymerase II has distinct roles in configuring nascent and steady state chromatin" [https://www.embopress.org/doi/full/10.1038/s44319-024-00085-x](https://www.embopress.org/doi/full/10.1038/s44319-024-00085-x).

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
