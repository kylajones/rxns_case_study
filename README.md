# rxns_case_study
KJ additions to separable hybrid models reaction engineering case study. Please run all notebooks in order presented below. Depdendencies and generated figures & sub-directories are listed in order of appearance in notebook.

# Notebook 1:
reaction_genExpt_v3_kj.ipynb

## Dependencies:
* files: 
  * reaction_models.py
  * lhs_20.csv
* sub-directories & files: 
  * training_data_small_noise_full_physics_2/
    *  
  * full_physics_kinetics/

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

# Notebook 3:
plateau_analysis.ipynb

## Dependencies (in order of appearance):
* files:
  * reaction_models.py
  * plateau_analysis_fxns.py
  * fitted_parameters_simple_physics_master.csv
 * sub-directories:
   * plateau_analysis_figs/
     * sp_fp_obj_plots/
       * case_i/, i = 1,...,10
     * sp_fp_plots/
       * case_i/, i = 1,...,10

## Generates:
* figures in:
  * plateau_analysis_figs/
    * sp_fp_obj_plots/
      * case_i/
        * fig_300_case_i.png
        * fig_350_case_i.png
        * fig_400_case_i.png
        * fig_450_case_i.png
        * fig_500_case_i.png 
    * sp_fp_plots/
      * case_i/
        * fig_300_case_i.png
        * fig_350_case_i.png
        * fig_400_case_i.png
        * fig_450_case_i.png
        * fig_500_case_i.png 
   * obj_vs_t_all_T_case_no_i_FP.png
   * obj_vs_t_all_T_case_no_i_SP.png
