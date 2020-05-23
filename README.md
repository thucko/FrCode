# FrCode
Analysis code for FrPNC experiments.

Spectra_fitter.py is a spectra fitting program. Currently it has two fitting functions, Lorentzain and Voigt. The current purpose of the program is to fit spectra from the E1 Stark induced transition. The Spectra_fitter_linux.py will be for linux distros until jax becomes available on windows.

Stark_shift_fitting.py reads in the Stark shift data file produced by Spectra_fitter_.py, and plot the frequency shift versus the electric field squared. It uses iminiuit to fit the data to a straight line.

peak_plotter.py is used for reading in the individual peak positions (and reduced chi squares (RCS)) from Spectra_fitter_.py and plotting them as a fucntion of time. It gives and y vs x plot and a histogram for the peak positions (RCS). Good for looking at the scatter in the peak positions during data acqusition.
  
Within the ~/Cavity is code for optical cavity parameters and light coupling calculations.
