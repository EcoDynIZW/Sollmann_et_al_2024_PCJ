# Data and code for "Accounting for observation biases associated with counts of young when estimating fecundity: case study on the arboreal-nesting red kite (Milvus milvus)", Sollmann, R., Adenot, N., Spakovszky, P., Windt, J., Mattsson, B.J. (submitted), PCJ.
--- 

All code necessary to implement the estimation simulations and population projections described in the above references manuscript

DOI repository: TBD2

DOI submitted pre-print: TBD2

## Description of the data and file structure

<ins> **This repository contains the following scripts (folder R)** </ins>:

*Parameter estimation simulations.R*: R code to generate imperfect and perfect counts (of nestlings) and fit the models described under Q1-Q3 in the above referenced manuscript, estimating average number of nestlings while accounting for false negative and false positive observation error

*Simulation summaries and plots.R*: R code to summarize simulation output and create plots

*Population projection.R*: R code to implement population projection described in the above referenced manuscript and produce plots

<ins> **It further contains the following model code files (folder R)** </ins>:

*Nimble model code.R*: Nimble model code for the models described in the manuscript and used in the estimation simulation

<ins> **It further contains the following data files (folder data-raw)** </ins>:

*datmat.rds*: Matrix with empirical observations of nestling counts, rows are true counts (climb counts), columns corresponding uncertain ground counts; values in cells correspond to number of nests with particular combination of true and uncertain count.

*Pnest.scenarios.rds*: Matrix with values of pi (probability a nest contains 1-4 nestlings); true input values, expected observed values based on uncertain (ground) counts only; estimated values from all scenarios investigated under Q3 (see main text for details)


## Sharing/Access information

NA


## Code/Software

Fitting models requires the R package Nimble. Please see publication for version information and R scripts for additional packages used in processing data. 
