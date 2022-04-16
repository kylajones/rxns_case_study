#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import os
import reaction_models as rxn
import pandas as pd

def get_case_params(case_no, sp_params):
    
    # load simple physics parameters
    θ_sp = sp_params.loc[sp_params['case'] == case_no] # select parameters for case
    θ_sp = θ_sp.drop(columns = 'case')
    θ_sp = θ_sp.to_numpy()
    θ_sp = θ_sp.reshape([4,1])
    return θ_sp


def spfp_plot_props(C_fp, C_sp, t, CA0, T, case_no, save):
    
    CA_fp, CB_fp, CC_fp = C_fp
    CA_sp, CB_sp, CC_sp = C_sp
    
    plt.figure()
    
    w_in = 3
    h_in = 4
    dpi_fig = 1200
    
    fig, ax = plt.subplots(figsize = (w_in, h_in))
    
    # plot full physics model
    plt.plot(t, CA_fp, 'b-', linewidth = 1, label = r'$C_A$')
    plt.plot(t, CB_fp, 'g-', linewidth = 1, label = r'$C_B$')
    plt.plot(t, CC_fp, 'r-', linewidth = 1, label = r'$C_C$')
    
    # plot simple physics model
    plt.plot(t, CA_sp, 'b--', linewidth = 1)
    plt.plot(t, CB_sp, 'g--', linewidth = 1)
    plt.plot(t, CC_sp, 'r--', linewidth = 1)
    
    # for tricking the legend
    # plt.plot(t[0], CA_sp[0], 'k-', label = r'FP')
    # plt.plot(t[0], CA_sp[0], 'k--', label = r'SP')
    
    # title, axis labels, legend
    plt_tit = r'$C_{A0}=$'+str(CA0)+r' [M], T='+str(T)[:3]+r' [K]'
    plt.title(plt_tit, fontsize = 14)
    plt.xlabel(r't [hr]', fontsize = 14)
    plt.ylabel(r'C [M]', fontsize = 14)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize = 12)
   
    # additional properties
    rxn.set_aspect_ratio(ax, 1.0)
    plt.xticks(np.arange(0.0,1.1,0.2))
    plt.xlim(0, 1)
    plt.ylim(0, 2)
    plt.yticks(np.arange(0.0, 2.5, 0.5))
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    if save == 'TRUE':
        cwd = os.getcwd()
        path_to_fig = os.path.join(cwd, 'plateau_analysis_figs', 'sp_fp_plots', 'case_'+str(case_no))
        file_name = 'fig_'+str(T)[:3]+'_case_'+str(case_no)+'.png'
        fig.savefig(os.path.join(path_to_fig, file_name), bbox_inches='tight', dpi= 1200)

def objective_function(CA, CB, CC, time, Temp, obj_coeff, print_level=2, temp_ref=500):
    '''
    inputs:
    CA: predicted concentration of species A; 2D numpy array or scalar
    CB: predicted concentration of species B; 2D numpy array or scalar
    CB: predicted concentration of species B; 2D numpy array or scalar
    time: prediction times; 2D numpy array or scalar
    Temp: prediction temperatures; 2D numpy array or scalar
    obj_coeff: coefficients or weights to use in objective function calculation 1x5 list / array
    [CA_coeff,CB_coeff,CC_coeff,temp_coeff,time_coeff]
    model_name: model which was used to generate predictions, to annotate console outputs
    print_level: verbosity of console outputs
    temp_ref: Reference value to normalize temperature in objective function. default is 500 K.
    returns:
    obj: objective values, 2D numpy array or scalar
    '''

    A_obj_coeff = obj_coeff[0]
    B_obj_coeff = obj_coeff[1]
    C_obj_coeff = obj_coeff[2]
    temp_obj_coeff = obj_coeff[3]
    time_obj_coeff = obj_coeff[4]

    # provide a warning in case any temperature exceeds the reference
    # temperature used to normalize the objective function
    if temp_ref < np.amax(Temp):
        print('\nWARNING: Reference temperature (temp_ref={0} K) to normalize objective function is less than maximum temperature in mesh grid. Consider increasing temp_ref for accurate results'.format(temp_ref))
    # End warning

    obj = A_obj_coeff*CA + B_obj_coeff*CB + C_obj_coeff*CC + temp_obj_coeff*Temp/temp_ref + time_obj_coeff*time

    return obj

# end objective_function 


# In[ ]:

def plot_obj_fxn(t, obj_sp, obj_fp, CA0, T, case_no, save):
    
    w_in = 3
    h_in = 4
    dpi_fig = 1200
    
    fig, ax = plt.subplots(figsize = (w_in, h_in))
    
    # plot full and simple physics objectives
    plt.plot(t, obj_sp, 'k--', label = r'SP')
    plt.plot(t, obj_fp, 'k-', label = r'FP')
    
    # title, axes, legend
    plt_tit = r'$C_{A0}=$'+str(CA0)+r' [M], T='+str(T)[:3]+r' [K]'
    plt.title(plt_tit, fontsize = 14)
    plt.xlabel(r't [hr]', fontsize = 14)
    plt.ylabel(r'obj [M]', fontsize = 14)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    # additional properties
    plt.xlim(0, 1.0)
    plt.ylim(-1.0, 2.0)
    plt.xticks(np.arange(0.0, 1.1, 0.2))
    plt.yticks(np.arange(-1.0, 2.5, 0.5))
    ax.tick_params(axis='both', which='major', labelsize = 12)
    
    rxn.set_aspect_ratio(ax, 1.0)
    
    if save == 'TRUE':
        cwd = os.getcwd()
        path_to_fig = os.path.join(cwd, 'plateau_analysis_figs', 'sp_fp_obj_plots', 'case_'+str(case_no))
        file_name = 'fig_'+str(T)[:3]+'_case_'+str(case_no)+'.png'
        fig.savefig(os.path.join(path_to_fig, file_name), bbox_inches='tight', dpi= 1200)

def plot_obj_fxn_temp(t, obj, CA0, T, case_no, save, model_type):
    
    w_in = 3
    h_in = 4
    dpi_fig = 1200
    
    fig, ax = plt.subplots(figsize = (w_in, h_in))
    
    # plot objective
    for j in np.arange(0, len(T)):
        if model_type == 'SP':
            plt.plot(t, obj[:,j], '--', label = r'T = ' + str(T[j])[:3] + r' K')
        else:
            plt.plot(t, obj[:,j], '-', label = r'T = ' + str(T[j])[:3] + r' K')
    
    # title, axes, legend
    plt_tit = r'Case '+str(case_no)+r': $C_{A0}=$'+str(CA0)+r' [M]'
    plt.title(plt_tit, fontsize = 14)
    plt.xlabel(r't [hr]', fontsize = 14)
    plt.ylabel(r'obj [M]', fontsize = 14)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize = 12)
    
    # additional properties
    plt.xlim(0, 1.0)
    plt.ylim(-1.0, 2.0)
    plt.xticks(np.arange(0.0, 1.1, 0.2))
    plt.yticks(np.arange(-1.0, 2.5, 0.5))
    ax.tick_params(axis='both', which='major', labelsize = 12)
    rxn.set_aspect_ratio(ax, 1.0)
    
    if save == 'TRUE':
        cwd = os.getcwd()
        path_to_fig = os.path.join(cwd, 'plateau_analysis_figs')
        file_name = 'obj_vs_t_all_T_case_no_'+str(case_no)+'_'+str(model_type)+'.png'
        fig.savefig(os.path.join(path_to_fig, file_name), bbox_inches='tight', dpi= 1200)
    