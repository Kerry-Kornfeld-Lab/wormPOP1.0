wormPOP
	The agent-based model wormPOP simulates the population dynamics of the nematode worm C. elegans cultured in a simple laboratory ecosystem. The laboratory ecosystem and wormPOP simulation include only two species – worms are the agents, and E. coli bacteria are the food source. wormPOP is designed to facilitate the investigation of how environmental factors (bacterial food and culling) and the life history traits of individuals (growth curves, fertility curves, mortality rates) impact population dynamics. wormPOP is designed to provide a detailed examination of individual worms based on explicit laboratory measurements of these behaviors, and aggregates the behaviors to predict population dynamics and demography. Specifically, wormPOP (1) tracks the life history of every individual worm (2) outputs data for the average behavior of individuals, and (3) combines individual behaviors to generate the emergent properties of population dynamics. 

This package includes the following files:
WormEnergy.dpr -- the main Delphi project file with almost all the Pascal code
WormEnergy.exe – Compiled, executable file for Windows
SysUtils.pas -- part of the Delphi7 package.
DosStuff.pas -- utilities for formatting text input and output
POISSONU.pas -- Generates random variables with a Poisson distribution
WormNRG.inp -- a text file in which input parameters are specified
summary_input_output_files.xlsx – summary and cross reference of inputs used for all figures, key for output files

System Requirements
	This PASCAL source is written to be compiled by Borland Delphi vx 7 (2002). With trivial modifications, I believe it can be made to run under Lazarus, and with a few more tweaks under Free Pascal. The main file is the DPr file, and it "uses" three PAS files, two of which are creations of Josh Mitteldorf, and the third, SysUtils.pas, is distributed with the Delphi Pascal compiler.
	It runs under Windows vx 7+ but operates in a DOS window, with only a text interface.
	Input is in a text file called WormNRG.inp, which is self-documenting.
	Output is in the form of CSV files.

Installation guide
	Download WormNRG.inp, WormEnergy.exe, on a PC with Windows vx 7+ in the same folder. Installation time should be under a minute. 

Instructions to use
	Specify the input parameter in Worm.NRG, save, and start WormEnergy.exe. The duration of the program is dependent on the chosen number of time steps- usually a few minutes. The outputting CSV files are saved in the same folder. For a new simulation, the existing files need to be removed or the program run in a new folder. 

Instructions to run on data
	Use the inputs in the provided example Worm.NRG file and run Worm.Energy.exe. It should take only a few minutes. The output file POPDYN shows the summary of the population at each time step (see output_key_POPDYN in summary_input_outout_files), the output file WRAtes reports all worm rates (wr- see Figure 1G), and the output file IndividualWorms shows the summary of single worms- but is empty and not used at this point. The CVS files POPDYN and WRates can then be processed in R, Excel, or different software.

Reproduction instructions
	In order to reproduce the data used in manuscript “A population dynamics tipping point for aging as a cause of adult death”- use the inputs listed in “input_files_all_figures” in the excel file “summary_input_output_files”.

Questions and licensing:
	This code was written by Josh Mitteldorf, and he will answer your questions if you email him at <aging.advice@gmail.com>. The executable and source code are released for non-commercial usage with attribution to the author. 

MIT Open Source License:
Copyright (c) 2021 by Josh Mitteldorf
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
