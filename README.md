# Pharmacokinetic-parameters-calculation
The calculation of the pharmacokinetic (PK) parameters using the concentration-time data.

# Description
This is a step by step guide of how PK parameters are computed.
This script can be used to validate the computation of PK parameters from the clinical software such as Phoenix WinNonlin or other programs used for this purpose

### What are the main PK parameters?
Cmax - the maximum plasma concentration or peak exposure, and the time to maximum plasma concentration
Tmax - time to achieve Cmax
Kel - elimination constant (how quickly the drug is eliminated from the bloodstream)
T1/2 - elimination half-life
AUC0-t - area under the curve of a plasma concentration versus time profile
AUClast - area under the curve from the time of dosing to the time of the last measurable (positive) concentration
AUCinf - area under the curve from time of dosing extrapolated to infinity, based on the last observed concentration
and others...


# Example of the calculation
The script takes the dataframe (excel sheet) and computes the sample size, power and pooledCV for the desired study
based on Cmax parameter (usually Cmax has the highest variability than AUCt).
Here is the diagram of how the calculation is going on:

1. Create an excel sheet with the input data of a single or multiple bioequivalence pivotal trials.

<img width="1280" height="720" alt="Presentation2" src="https://github.com/user-attachments/assets/19080d9a-3778-4246-b7bd-0c978535a61e" />

2. Run the script.
   
<img width="1280" height="720" alt="Presentation1" src="https://github.com/user-attachments/assets/9cf86d18-e531-4210-917c-10ababe68626" />

# Installation instructions
Please go to file named sample size pooled.R and click in the upper right corner the icon copy to copy the code, insert it, for example, into R studio,
specify the path for your excel file and run it.

# Found a bag?
If you found a bag or have an idea of how to improve this project, please contact me on [here](https://www.linkedin.com/in/vlia/) 

# What was learned from this project?
Sample size can be calculated by the central and the non central t approximation. 
The calculation based on the large sample approximation should be avoided because this underpowers the study.

# Ackowledgments
Some parts of the sample size calculation functions are adapted from the Helmut Schutz article regarding the comparison of the sample size calculation using
different approaches: large sample approximation, central t approzimation and noncentral t approximation
Please check this article [here] (https://bebac.at/articles/Sample-Size-Estimation-for-Equivalence-Studies-in-a-Parallel-Design.phtml)

Extended by: Valery Liamtsau
