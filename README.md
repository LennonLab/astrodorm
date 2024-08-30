# astrodorm
===================

A simple individual-based model to explore these ideas where an autocatalytic molecule replicates and decays while stochastically transitioning between active and inactive states in the presence of catalyst. This repository contains open-source code, data, and figures generated. 

### Repo Contents

* **figures:** - Figures generated from the code
  * Fig4a.png - A png file describing the how the proportion of active molecules that vary in their stochastic transition to dormancy vary as the abundance of catalyst
  * Fig4b.png - A png file describing the log_10(Abundance) of the autocatalytic molecule as its decay rate in dormancy varies

* **code:** - Analytical and figure generating code
  * NewMod_20240723_1413_AvgPop.R - Code to generate the text file to examine the decay of an autocatalytic molecule with dormancy
  * NewMod_20240723_1413_RatioFig.R - Code to generate a figure examining the ratio of active molecules to total molecules as stochastic transitions to dormancy and the number of catalysts vary
  * PopFig_20240727_2309.R - Code to generate a figure of the abundance of the autocatalytic molecule as its decay rate in dormancy varies
  * PopFig_20240809.Rmd - Code to generate a figure of the log_10(abundance) of the autocatalytic molecule as its decay rate in dormancy varies
  * Ratio_20240809.Rmd - Code to generate a figure examining the ratio of active molecules to total molecules as stochastic transitions to dormancy and the number of catalysts vary
