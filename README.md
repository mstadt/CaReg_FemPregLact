# Calcium regulation
 Code for mathematical models of calcium homeostasis in a male rat, female (virgin) rat, pregnant rat, and lactating rat. 
 This code was used to generate the results in ["Mathematical modeling of calcium homeostasis in female rats: An analysis of sex differences and maternal adaptations" by Melissa Stadt & Anita Layton](https://www.sciencedirect.com/science/article/pii/S0022519323001807).
 Copies of the figures are found under "FinalFigures/"
 
## Files for making figures in manuscript

**Figures 3.1 and 3.2** are made using **makefig_SS.m**. This file uses steady state results that can be created using **compute_1SS.m**

**Figure 3.3** is made using **makefig_localsens_all.m**. This file uses local sensitivity results creased from **compute_localsensitivity.m**

**Figure 3.4** is made using **makefigs_male2fem_fem2male.m**. This file uses results from **compute_male2fem_fem2male_sensitivity.m**

**Figure 3.5** is made using **makefigs_vir2preg_vir2lact_sensitivity.m**. This file uses results from **compute_vir2preg_vir2lact_all.m**

**Figure 3.6** is made using **makefigs_preg2fem_lact2fem_sensitivity.m**. This file uses results from **compute_preg2fem_lact2fem_all.m**

**Figure 3.7** is made using **makefigs_preg2lact_lact2preg_sensitivity.m**. This file uses results from **compute_preg2lact_lact2preg.m**

**Figures 3.8 & 3.9** is made using **makefigs_D3deficient.m**. This file uses results from **driver_D3deficiency.m**

**Figure 3.10** is made using **SS_malepreglact.m**. This file uses results from **compare_male2preglact.m**

**Figure 3.11** is made using **SS_malepreglact_regcalcitriol.m**. This file uses results from **compare_male2preglact.m**

## Key files

**params_male.txt, params_female.txt, params_pregnancy.txt, params_lactation.txt** list the parameter values for the male rat, female rat, pregnant rat, and lactating rat models

**calcium_mod.m** contains the model equations
