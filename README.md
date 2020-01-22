# FrCode
Analysis code for FrPNC experiments.

FrPNC_Stark_v1.py is a spectra fitting program. Currently it has two fitting functions, Lorentzain and Voigt. The current purpose of the program is to fit spectra from the E1 Stark induced transition.

Next will be to add plotting of the frequency shifts and output the final result. The will require:
  1) Saving off peak postions 
  2) conversion of applied voltage to electric field strength
  3) fiiting of either a quadratic or linear function
  
The spectra fitting portion of the code can be used for other experiments. I plan to adapt the spectra fitting for M1 and improve its speed (if possible).
