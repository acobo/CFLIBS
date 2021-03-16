# CFLIBS
Calibration-free LIBS algorithm in Matlab

This is a Matlab implementation of the CF-LIBS algorithm originally proposed by Ciucci et al: (CIUCCI, A., et al. New procedure for quantitative elemental analysis by laser-induced plasma spectroscopy. Applied spectroscopy, 1999, vol. 53, no 8, p. 960-964).

This implementation tries to get a universal algorithm (as much as posible!) with the following features:

1. it is a long single-file code with almost everything needed (database of emission lines and spectroscopic parameters, expected composition of the sample, different processing steps and algorithms... The external dependences are just a few, for example, the captured spectrum of a calibration white-light source to know the spectral response of the experimental setup.
2. At this moment, it cannot be executed because we are including all the dependences and there is no sample data to play with. yet. sorry.
3. Uses a user-specified a-priori selection of lines for each element expected in the sample. Some of the emission lines could (and should) be automatically discarded by the algorithm further in the processing chain.
4. The specific processing steps, parameters and algorithms are specified using a set of configuration variables.
5. under construction...

