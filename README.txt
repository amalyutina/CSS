R scripts:

1. calculation_of_css.Rmd
 The script uses single_drug_data.RData and combo_drug_data.RData to calculate drug sensitivity scores (DSS) for single drugs and combination sensitivity scores (CSS) for drug combinations. The resulting CSS are saved to the css_data.RData file. 

2. predicion_of_css.Rmd
 The script predicts the CSS for drug combinations from css_data.RData using drug target information from fingerprints.RData, primary_targets.RData and primary_plus_sea_targets.RData, and machine learning methods (Elastic Net, Random Forests, Support Vector Machines). R2, RMSE, MAE and correlation values are calculated for a selected cell line and machine learning method.

RData files:

1. single_drug_data.RData
 A file that can be used in  calculation_of_css.Rmd to calculate DSS for single drug experiments. The file has to contain such columns as cellLine, drug, dose and viability. Column batchID is optional and it is possible to have multiple viability columns (replicates), but their names should start with ‘viability’.

2. combo_drug_data.RData
A file that can be used in  calculation_of_css.Rmd to calculate CSS for combination drug experiments. The file has to contain such columns as cellLine, drugA, drugAdose, drugB, drugBdose and viability. Column batchID is optional and it is possible to have multiple viability columns (replicates), but their names should start with ‘viability’.

3. fingerprints.RData
A file that can be used in predicion_of_css.Rmd to build a matrix of features. It contains binary fingerprint data mapped to the drugs used in the experiment. It has to contain column  with drug names named ‘drug’.

4. primary_targets.RData
A file that can be used in predicion_of_css.Rmd to build a matrix of features. It contains binary primary drug target data mapped to the drugs used in the experiment. It has to contain column  with drug names named ‘drug’.

5. primary_plus_sea_targets.RData
A file that can be used in predicion_of_css.Rmd to build a matrix of features. It contains binary primary drug target data and drug targets that were identified using SEA method mapped to the drugs used in the experiment. It has to contain column  with drug names named ‘drug’.

6. css_data.RData
A file that can be used in predicion_of_css.Rmd. The file provides CSS data that were calculated using calculation_of_css.Rmd.  It has to contain columns  with both drug names named ‘drug1’ and ‘drug2’, cell lines specified in column ’cellLine’ and CSS column names ‘CSS’.