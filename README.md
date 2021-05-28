# TCS_analysis_CLAS12

This repository contains some of the scripts used for the TCS analysis. There are in there "raw" form, with dirty formatting.

analysisTCSn1CheckSystematicsWithNewaccalgo1Dacc : main analysis code
run as clas12root -l analysisTCSn1CheckSystematicsWithNewaccalgo1DaccSimu.C DATA_FILE Acceptance_FILE SEED
The seed is unused, insert 1, for example

pippim : used to derive the efficiency and central momentum correction

Systematic.pl/NiceSytematic.pl/Plot.C/PlotNoSyst.C : perl scripts running the systematic combination and doing the final plotting
