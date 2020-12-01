# FrCode
Analysis code for FrPNC experiments.

spectra_fitter.py is a spectra fitting program. You can fit to two of the most common fucntions used in spectroscopy, Lorentzain and Voigt. 


Stark_shift_fitting.py reads in the Stark shift data file produced by Spectra_fitter_.py, and plot the frequency shift versus the electric field squared. It uses iminiuit to fit the data to a straight line.

peak_plotter.py is used for reading in the individual peak positions (and reduced chi squares (RCS)) from Spectra_fitter_.py and plotting them as a fucntion of time. It gives and y vs x plot and a histogram for the peak positions (RCS). Good for looking at the scatter in the peak positions during data acqusition.
  
Within the ~/Cavity is code for optical cavity parameters and light coupling calculations.

# Note:
The plotting style used is a custom style and will not work on other machines.
