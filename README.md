# TCS_analysis_CLAS12

This repository contains some of the scripts used for the TCS analysis. There are in there "raw" form, with dirty formatting.

___
analysisTCSn1CheckSystematicsWithNewaccalgo1Dacc : main analysis code
run as 

clas12root -l analysisTCSn1CheckSystematicsWithNewaccalgo1DaccSimu.C DATA_FILE Acceptance_FILE SEED

The seed is unused, insert 1, for example
___
pippim : used to derive the efficiency and central momentum correction
___
Systematic.pl/NiceSytematic.pl/Plot.C/PlotNoSyst.C : perl scripts running the systematic combination and doing the final plotting

___
Extensive documentation on this analysis is found in:
https://www.jlab.org/Hall-B/shifts/admin/paper_reviews/2021/AnalysisNoteV3_6-1461061-2021-04-09-v8.pdf
