# drug_treatment_thesis

This repository is used to summarise the source code I wrote for my thesis (to be completed in Autumn 2022): Multiparametric analysis of drugs on YAP localization. The project aims to used the observational data to analyse the effect of 15 drugs on localisation of YAP (Yes-associated protein) a protein that directs gene expression. The project is under supervision of Dr Sarah Flippi and is a extended work to Dr J.Sero's article: Multiparametric Analysis of Cell Shape Demonstrates that β-PIX Directly Couples YAP Activation to Extracellular Matrix Adhesion.

## General information
In the analysis, Random Forest is the main model will be used. Both regression and classification RF model will be applied to model the 1) YAP localisation corresponds to other cell features and 2) classify the treatment/control data by cell features including YAP ratio.


## current avaliable files

`preprocessing.R`: deal with data pre-processing
`T-learner.R`: Inspired by article Metalearners for estimating heterogeneous treatment effects using machine learning (Künzel et al., 2019), an function for analysing the causal effect of treatment to response variable.
`hparan tuning.Rmd`: file used for tune the hyper-parameter for Random Forest models.
`functions.R`: support functions including parallel operation of hyper-parameter tuning.

