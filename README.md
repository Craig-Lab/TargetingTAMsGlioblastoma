# TargetingTAMsGlioblastoma

### Retrieving virtual clinical trial datasets
Data sets VCT_none.mat , VCT_SOC.mat and VCT_SOC_ICI.mat are found on dataverse: https://dataverse.harvard.edu/privateurl.xhtml?token=78088a14-bd0e-44a8-b5a6-fb8a35411a67 . They need to be included in the VCT folder. Do not change their name.

### Reproducing figures
The file Panels.m contains all the code necessary to reproduce the figures in the article. Do not change the hierarchy of the folders if you want the code to work as it is. 

### Running simulations
There are 3 files necessary to run simulations:
- set_patients.m
- set_treatment.m
- solver_gbm.m

Run the files one after another, in that order. 

set_patients.m sets the values of the main model parameters (non-treatment related) and of those needed to run the solver if you want to simulate a tumour growing from initiation without treatment (omit set_treatment.m then).

set_treatment.m sets the treatment schedule. 

Before running the solver, you can adjust the following boolean parameters to decide what treatments to simulate (0=not administered, 1=administered):

- p.treatmentBoolean: resection
- p.treatmentTMZ: chemotherapy
- p.treatmentRT: radiotherapy
- p.treatmentICI: nivolumab
- p.treatmentMac: TAM-targeting strategies

solver_gbm.m contains the ODE solver. The function fullmodel(t,y,p) is called when no ICI is administered, whereas the function fullmodelICIadmin(t,y,p,initialTime) is called after an ICI dose is administered since the ICI ODE is dependant on whether a dose is given or not.

### Supplementary questions
For any questions regarding the code, please write to blanche.mongeon@umontreal.ca .


