=============================================================================
geostats2D: A MatLab library for geostatistical analysis of 2D datasets
=============================================================================

This repository contains codes developed in MATLAB programming environment for the purposes of my Diploma thesis in School of Mineral Resources Engineering of Technical University of Crete. These codes are relevant to  Geostatistics and more specifically to variographic analysis and implementation of Ordinary Kriging estimation in 2D datasets.
=============================================================================

This repository contains the following files:

randomfield.m : constructs random field with specified covariance model
detrendv.m : detrends 2D data with various trend models using linear regression
ccvfun.m : calculates the correlation, covariance and variogram for given model and distances
expvar.m : calculates the experimental (directional or omnidirectional) variogram
variogramfit.m : fits a desired model to the given experimental variogram
aniso_dvf.m : estimates the anisotropy parameters of a 2D dataset by means of directional variogram fitting
crossval.m : calculates the cross-validation scores of Ordinary Kriging performance of a given model
ordkrig.m : implements ordinary kriging with a specified model for estimating missing values of a field 
example1.m, example2.m : examples of the above codes usage
 
For more information see help block in functions
=============================================================================
It is strongly recommended not to assume that these codes are entirely correct, and in no case I claim them as such. I am publishing these codes only for giving to the interested a basis for variography and ordinary kriging implementation and in the hope that it will be 
useful. You can use and modify these codes. Comments and questions are welcome.