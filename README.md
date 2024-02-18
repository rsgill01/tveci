# tveci
Code for fitting Threshold Error Correction model with Data Imputation

This repository includes the following five files.

1. **Vignette.pdf** demonstrates how to use the functions `c.tveci` and `r.tveci` in the R statistical computing environment to estimates the parameters of the threshold error correction model for data on daily minimum gas prices in Perth based on prices on the previous seven days and the crude oil price on the previous day.

2. **PerthMin2015.dat** is a tab-delimited file with the data used the example described in the vignette.

3. **Rcode.R** is a text file with source code written in the R language.

4. **tvec.c** gives some C functions called by the R code in **Rcode.R**.

5. **tvec.dll** is a dynamic-link library with the compiled C code from **tvec.c**.
