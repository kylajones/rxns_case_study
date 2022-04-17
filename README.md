# rxns_case_study
KJ modifications and additions to separable hybrid models reaction engineering case study. Please run all notebooks in order presented below. Note that indices i and j represent case study numbers (i = 1,...,10) and experiment numbers (j = 1,...,19), respectively.

# Notebook 1:
reaction_genExpt_v3_kj.ipynb

## Dependencies:
* files: 
  * reaction_models.py
  * lhs_20.csv
* sub-directories & files (if running with EE data, otherwise these files are created for a new user): 
  * training_data_small_noise_full_physics_2/
    *  A_true.csv
    *  E_true.csv
    *  rxn_exp_{j}.csv
    *  rxn_exp_{j}_nonoise.csv
  * full_physics_kinetics/

## Generates:
* file: fig-rxn_lhs.png
* files in sub-directory:
  * full_physics_kinetics/
    * fig-exp_{j}.png

# Notebook 2:
sp_fitted_params.ipynb

## Dependencies:
* files: 
  * reaction_models.py
  * experiment_list.csv
* sub-directories & files: 
  * regressed_params_sp/
    * fitted_parameters_simple_physics_{i}.csv
  * training_data_small_noise_full_physics_2/
    *  A_true.csv
    *  E_true.csv
    *  rxn_exp_{j}.csv
    *  rxn_exp_{j}_nonoise.csv
* sub-directories:
  * simple_physics_plots/
    * case_{i}/

## Generates:
* files: 
  * in simple_physics_plots/case_{i}/
    * case_i_exp_{j}.png
  * A1_vs_case.png
  * A2_vs_case.png
  * E1_vs_case.png
  * E2_vs_case.png

# Notebook 3:
plateau_analysis.ipynb

## Dependencies:
* files:
  * reaction_models.py
  * plateau_analysis_fxns.py
  * fitted_parameters_simple_physics_master.csv
 * sub-directories:
   * plateau_analysis_figs/
     * sp_fp_obj_plots/
       * case_{i}/
     * sp_fp_plots/
       * case_{i}/

## Generates:
* figures in:
  * plateau_analysis_figs/
    * sp_fp_obj_plots/
      * case_{i}/
        * fig_300_case_{i}.png
        * fig_350_case_{i}.png
        * fig_400_case_{i}.png
        * fig_450_case_{i}.png
        * fig_500_case_{i}.png 
    * sp_fp_plots/
      * case_{i}/
        * fig_300_case_{i}.png
        * fig_350_case_{i}.png
        * fig_400_case_{i}.png
        * fig_450_case_{i}.png
        * fig_500_case_i.png 
   * obj_vs_t_all_T_case_no_{i}_FP.png
   * obj_vs_t_all_T_case_no_{i}_SP.png
