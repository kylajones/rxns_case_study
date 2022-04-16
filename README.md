# rxns_case_study
Separable hybrid models for reaction engineering case study
add a description of all files


# Notebook 1:
reaction_genExpt_v3_kj.ipynb

## Dependencies (in order of appearance):
* file: reaction_models.py
* file: lhs_20.csv
* sub-directory & files (included in repo): training_data_small_noise_full_physics_2
* sub-directory: 'full_physics_kinetics'

## Generates:
* file: fig-rxn_lhs.png
* files in sub-directory : cwd+'full_physics_kinetics/fig-exp_{i}.png' i = 0,...,19

# Notebook 2:
sp_fitted_params.ipynb

## Dependencies (in order of appearance):
* file: reaction_models.py
* sub-directory & files (included in repo): 'regressed_params_sp'
* sub-directory & files (included in repo): 'training_data_small_noise_full_physics_2'
* file : experiment_list.csv
* sub-directory/sub-sub-directories : 'simple_physics_plots/case_i' i = 1,..,10

## Generates:
* files: 
  * 'case_i_exp_j.png', i= 1,...,10, j=0,...,19
  * 'A1_vs_case.png'
  * 'A2_vs_case.png'
  * 'E1_vs_case.png'
  * 'E2_vs_case.png'
