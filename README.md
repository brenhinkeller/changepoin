# changepoint

This repository contains a C implementation of 
transdimensional Bayesian changepoint analysis based on 
the approach of Gallagher et al. 2011 "Inference of 
abrupt changes in noisy geochemical records using trans-
dimensional changepoint models",
(doi.org/10.1016/j.epsl.2011.09.015)[https://doi.org/10.1016/j.epsl.2011.09.015)

The files gauss.h and pcg_variants.h contain code written
by others and previously released under open-source GNU 
and Apache licenses. All new code herein is released under 
a GNU GPL v2.0 licencse (see LICENSE)

## Usage

```
gcc -o changepoint changepoint.c
./changepoint 10 1000000 exampledata.csv > examplechangepoints.csv
```
