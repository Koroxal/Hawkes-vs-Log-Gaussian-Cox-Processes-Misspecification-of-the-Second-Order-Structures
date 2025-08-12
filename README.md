This repository contains the code associated with the paper:
"Misspecification issues between competitive spatio-temporal cluster point processes"
by the authors Alba Bernabeu, Claudio Fronterrè, and Jorge Mateu.

It provides generic code for the simulation and estimation of Hawkes processes and Log-Gaussian Cox Processes (LGCP).
For clarity and concision, only one case study simulation example is presented for each of Sections 3.1 and 3.2 of the paper.
For further information, please contact the authors (bernabeu@uji.es).

The folder sim_hawkes_estimate_lgcp contains the code associated with Section 3.1, representing a case of model misspecification where data are simulated from a Hawkes process, but second-order characteristics are estimated and analysed under the assumption that the data originate from a LGCP.

The folder sim_lgcp_estimate_hawkes corresponds to the code for Section 3.2, where data are simulated from a LGCP, but second-order characteristics are estimated and analysed under the assumption that the data come from a Hawkes process, representing the converse case of model misspecification.
