# INFERNO_eCAVIAR

This is a custom version of the colocalization method eCAVIAR (https://github.com/fhormoz/caviar), used for analyses with the INFERNO pipeline and related. 
This version:
* Outputs region-level colocalization probabilities, i.e. the probability that 0, 1, 2 variants in the region are colocalized.
The original only outputs SNP-level colocalization probabilities, and region-level probabilities for each of the individual studies.
* Implements some optimizations in the computation of Bayes factors, to improve runtime and numerical stability.

So far, this version only allows a maximum of 1 or 2 causal variants in a region.

## Installation
```
git clone https://github.com/matei-ionita/INFERNO_eCAVIAR
cd INFERNO_eCAVIAR
make
```

INFERNO_eCAVIAR depends on the Armadillo library for linear algebra computations; it is included in the repository for convenience.

## Usage
```
INFERNO_eCAVIAR -z <GWAS file> -z <eQTL file> -l <LD file> [options]
options:
  -z file containing a column of variant names and a column of variant summary statistics
  -l file containing LD matrix, tab-separated
  -c maximum number of causal variants allowed (default 2)
  -o name of output file
```
