This package runs in CONSOLE mode under Borland Delphi vx 7 (2002). With trivial modifications, I believe it can be made to run under Lazarus, and with a few more tweaks under Free Pascal.
It runs under Windows vx 7+ but operates in a DOS window, with only a text interface.
Input is in a text file called WormNRG.inp, which is self-documenting.
Output is in the form of CSV files.

This package includes the following files.

WormEnergy.dpr -- the main Delphi project file with almost all the Pascal code
DosStuff.pas -- utilities for formatting text input and output
POISSONU.pas -- Generatles random variables with a Poisson distribution
WormNRG.inp -- a text file in which input parameters are specified

This code was written by Josh Mitteldorf, and he will answer your questions if you email him at <aging.advice@gmail.com>.
It is released for non-commercial usage with attribution to the author, and permission to modify is granted, again for non-commercial uses.   