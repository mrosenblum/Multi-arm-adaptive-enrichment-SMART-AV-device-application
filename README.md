# Multi-arm-adaptive-enrichment-SMART-AV-device-application
Optimized Adaptive Enrichment Designs for Multi-Arm Trials: Learning which Subpopulations Benefit from Different Treatments, by Jon Arni Steingrimsson, Joshua Betz, Tianchen Qian, and Michael Rosenblum

The code in this directory runs the trial designs compared in the "Optimized Adaptive Enrichment Designs for Multi-Arm Trials: Learning which Subpopulations Benefit from Different Treatments" manuscript. To run the code, you need the following R packages:
mvtnorm,
plyr,
ggplot2,
gridExtra,
xtable
Once you've installed these and loaded them via the call:
library(mvtnorm,plyr,ggplot2,gridExtra,xtable)
you can run (i.e., source) the designOptimize2Treatments.R file, which will solve the 8 optimization problems from the paper (and will produce the 8 optimized designs). The time to run this is substantial.
