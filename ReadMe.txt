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
Supplemental_Material--information for wormPop1.0 design and use
Figure 2_lifespan_analysis --R code for Figure 2 M,N,P
Figure 4_egg_laying_1 -- R code for Figure 4F
Figure 4_egg_laying_2 -- R code for Figure 4F


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
	In order to reproduce the data used in manuscript “A laboratory and simulation platform to integrate individual life history traits and population dynamics”- use the inputs listed in “input_files_all_figures” in the excel file “summary_input_output_files”. Data for Figure 2 M,N,P were analyzed using R code “Figure 2_lifespan_analysis” and data for Figure 4F were analyzed using R code “Figure 4_egg_laying_1” and  “Figure 4_egg_laying_2”.
Details for using “Figure 2_lifespan_analysis”: Open RStudio (Version 1.4.1717), set working directory, open file “Figure 2_lifespan_analysis”, create a folder graphs in your working directory, save lifespan data according instruction in the code comments under file name “lifespan”, adjust line 92 candidate names to analyze, adjust line 151 according number of groups/treatments, and run code.
Details for using “Figure 4_egg_laying_1” and “Figure 4_egg_laying_2”: Open RStudio (Version 1.4.1717), open file “Figure 4_egg_laying_1” and “Figure 4_egg_laying_2”, set working directory, change in code your file_location (line 5) “Figure 4_egg_laying_1”, run code, store number of analyzed individuals and total number of eggs in new file for all expriments, save file and add file location to “Figure 4_egg_laying_2” in line 11. Run “Figure 4_egg_laying_2”. 

Questions and licensing:
	This code was written by Josh Mitteldorf with the exception of Rcode, and he will answer your questions if you email him at <aging.advice@gmail.com>. The executable and source code are released for non-commercial usage with attribution to the author.

MIT Open Source License:
Copyright (c) 2021 by Josh Mitteldorf
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



