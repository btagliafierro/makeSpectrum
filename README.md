# makeSpectrum
This class can be used to create elastic response spectrum from a ground motion record in MATLAB.
The file main.m reports the template to use the class to generate the response spectrum of any given signal. The provided example uses the outdated record of the KOBE event 1995, which uses a quite coarse sampling frequency.

The number of steps that compose the spectrum can be adjusted directly in the class definition (no set functions). The provided version uses variable time steps, with smaller dt for period that are more meaningfull for buildings.
