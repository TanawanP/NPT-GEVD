# R Code for NPT_GEVD_main file --------------------------------------------------------------------
# Title: Analysis of Negative Power Transformation-GEV Distribution (NPT-GEVD) for Estimation of Minimum IAT Data
# Authors: Prahadchai T., Yoon S., and Busababodhin P.
# Version 1 (1 Dec 2024)
# --------------------------------------------------------------------------
# Short Description:
# There are a total of 8steps to obtain the parameter and return level (RL) estimates.
# Please follow steps 1)–8). Note that steps 3)–7) require setting the directory, providing minimum data, and obtaining the results.
# Step 8) involves saving the final output in Excel format and storing it in your folder directory.
# --------------------------------------------------------------------------

# Step 1) Packages required
# Step 2) Important functions: 
## Function1 - Calculation function of the RL for CT-GEVD using the MLE method and the standard error (SE) of RL using the delta method.
## Function2 - Calculation function of the RL for RT-GEVD using the MLE method and the SE of RL using the delta method.
## Function3 - Calculation function of the RL for NPT-GEVD using the MLE method and the SE of RL using the delta method.
## Function4 - Bootstrap sampling technique for refining the SE of parameter estimates in cases where the Fisher information matrix returns NA values.
## Function5 - Calculation function of the RL for CT-GEVD using the LM method and the SE of RL using the bootstrap sampling method.
## Function6 - Calculation function of the RL for RT-GEVD using the LM method and the SE of RL using the bootstrap sampling method.
## Function7 - Calculation function of the RL for NPT-GEVD using the LM method and the SE of RL using the bootstrap sampling method.
## Function8 - Optimization function based on the L-BFGS-B method by minimizing the Cramér-von Mises (CvM) statistical test with LM method
## Function9 - Optimization function based on the L-BFGS-B method by minimizing the CvM statistical test with MLE method
## Function10 - Hyperparameter estimation function using resampling with an added small noise technique.

# Step 3) Setting path
# Step 4) Read a list of files in a folder 
# Step 5) The main function estimates parameters and RL
# Step 6) This function is required for converting the output from the function 'cal.g' to adjust the results table into a dataframe format.
# Step 7) Loop for obtaining output, where i,…,n represent the station indices.
# Step 8) Save output as an excel format 
