# FrCode
Analysis code for FrPNC experiments.

FrPNC_Spectra_fitter_v5.py is a spectra fitting program. Currently it has two fitting functions, Lorentzain and Voigt. The current purpose of the program is to fit spectra from the E1 Stark induced transition.

FrPNC_Stark_shift.py reads in the Stark shift data file produced by FrPNC_Spectra_fitter_v5.py, and plot the frequency shift versus the electric field squared. It uses iminiuit to fit the data to a straight line.

peak_plotter.py is used for reading in the individual peak positions (and reduced chi squares (RCS)) from Fr_PNC_Spectra_fitter_v5.py and plotting them as a fucntion of time. It gives and y vs x plot and a histogram for the peak positions (RCS). Good for looking at the scatter in the peak positions during data acqusition.
  
