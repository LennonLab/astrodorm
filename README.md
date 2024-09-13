# astrodorm
===================

A simple individual-based model to explore these ideas where an autocatalytic molecule replicates and decays while stochastically transitioning between active and inactive states in the presence of catalyst. This repository contains open-source code, data, and figures generated. 

### Repo Contents

* **figures:** - Figures generated from the code
  * Fig4a.png - A png file describing the how the proportion of active molecules that vary in their stochastic transition to dormancy vary as the abundance of catalyst
  * Fig4b.png - A png file describing the log_10(Abundance) of the autocatalytic molecule as its decay rate in dormancy varies

* **code:** - Analytical and figure generating code
  * Ratio+Persistence_20240912.Rmd - The code used to create a model where an autocatalytic molecule can enter and exit dormancy and interact with a catalyst. Used also to generate figures. 
  * Ratio.csv - The csv file output output to examine how the proportion of active autocatytic molecules varies in the presence of a catalyst. Used to make Fig4a.png
  * AverageShielding.csv - The csv file output to examine how decay in dormancy influences persistance. Used to make Fig4b.png
