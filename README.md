# *sc*RBA: *S. cerevisiae* Resource Balance Analysis model
This repository provides resources and model files for the genome-scale resource balance analysis (RBA) model *sc*RBA that accompanies the following manuscript:

  - Hoang V. Dinh and Costas D. Maranas. "Evaluating proteome allocation of *Saccharomyces cerevisiae* phenotypes with resource balance analysis." Metab. Eng., 2023, 77:242-255. https://doi.org/10.1016/j.ymben.2023.04.009

which have been updated to follow the conventions in the following manuscript: 
  - Eric J. Mooney, Patrick F. Suthers, Wheaton L. Schroeder, Hoang V. Dinh, Xi Li, Yihui Shen, Tianxia Xiao, Catherine M. Call, Heide Baron, Arjuna M. Subramanian, Daniel R. Weilandt, Felix C. Keber, Martin Wühr, Joshua D. Rabinowitz, Costas D. Maranas. "Metabolic flux and resource balance in the oleaginous yeast Rhodotorula toruloides." Metab. Eng. https://doi.org/10.1016/j.ymben.2025.11.012.

Please see the manuscripts above for details on the model, and the repository https://github.com/maranasgroup/rtRBA.

For the 2023 version, formulation and software usage guide are available at [suppMat/scRBA_suppText1_2023-04-10.docx](suppMat/scRBA_suppText1_2023-04-10.docx) and [suppMat/scRBA_suppText2_2023-04-10.docx](suppMat/scRBA_suppText2_2023-04-10.docx)
Program and packages requirements for the provided software implementation: Python 3 (plus packages: cobra, pandas, numpy, matplotlib, jupyter, scikit-learn, cplex (compiled from installed CPLEX linear programming solver), and all automatically-installed associated dependencies to those already listed) and GAMS with ("built-in") Soplex solver. (Tested on Python 3.6 and GAMS 39.1.0. Python package installation could be done through "pip" and could take up to 1 hour)