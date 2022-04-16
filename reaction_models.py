### This script contains the mathematical models for the cannonball problem and associated helper functions
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate
import scipy.optimize as optimize

#---------------------------------------------------
#KJ Checkpoint
#matplot lib error: unable to connect to server
#---------------------------------------------------
# helper function for plotting
# modified from https://stackoverflow.com/questions/27579898/how-to-create-square-log-log-plots-in-matplotlib
def set_aspect_ratio(plot, aspect_ratio):
        x_min, x_max = plot.get_xlim()
        y_min, y_max = plot.get_ylim()
        return plot.set_aspect(aspect_ratio * ((abs(x_max - x_min)) / (abs(y_max - y_min))))

def kinetics(A, E, T):
    ''' Computes kinetics from Arrhenius equation

    Arguments:
        A: pre-exponential factor, [1 / hr]
        E: activation energy, [kJ / mol]
        T: temperature, [K]

    Returns:
        k: reaction rate coefficient, [1/hr] or [1/hr*L/mol]

    '''
    R = 8.31446261815324 # J / K / mole

    return A * np.exp(-E*1000/(R*T))
#-----------------------------------------------------
#KJ Checkpoints
#A = np.array([200,100,50])
#E = np.array([10, 20, 15])
#T = 473
#k = kinetics(A, E, T)
#print('Checkpoint: k is calculated:',k)
#-----------------------------------------------------

def full_physics_model_2(theta,t,CA0,T):
    '''Evaluate model B numerically (adapted from NSA Spring 2020 solutions)

    Arugments:
        t: time, [hour], Numpy array
        theta: fitted parameters, Numpy array
        CA0: initial concentration, [mol/L], scalar
        T: temperature, [K], scalar

    Returns:
        CA, CB, CC: Concentrations at times t, [mol/L], three scalars or numpy arrays
    '''

    ### BEGIN SOLUTION

    # units: [1/hr]
    k1 = kinetics(theta[0], theta[3], T)

    # units: [1/hr * L/mol] --------------> recheck these units
    k2 = kinetics(theta[1], theta[4], T)

    # units: [1/h * l/mol]
    k3 = kinetics(theta[2], theta[5], T)

    def rhs(t_,y_):
        '''RHS of differential equation'''
        dy_ = np.zeros(3)

        # units: mol/L
        rA = -k1*y_[0] + k3*y_[1]
        rB = k1*y_[0] - k2*(y_[1]**3) - k3*y_[1]
        #rB = k1*y_[0] - k2
        rC = k2*(y_[1]**3)

        dy_[0] = rA
        dy_[1] = rB
        dy_[2] = rC

        return dy_

    # Integrate with solve_IVP
    tspan = [min(t), max(t)]
    y0 = [CA0, 0.0, 0.0]
    results = integrate.solve_ivp(rhs, tspan, y0, method="BDF",t_eval=t)

    # Extract solution
    CA = results.y[0,:]
    CB = results.y[1,:]
    CC = results.y[2,:]

    return CA, CB, CC

    ### END SOLUTION

#----------------------------------------------------------------
#KJ checkpoint
#theta = np.hstack([A,E])
#t = np.linspace(0,5)
#CA0 = 2
#C = full_physics_model_2(theta, t, CA0, T)
#print('full physics model 2 checkpoint: C =', C)
#----------------------------------------------------------------

# nonlinear parameter estimation with full physics model 2
def regression_func_full_physics_model_2(theta, data):
    '''
    Function to define regression function for least-squares fitting
    Arguments:
        theta: parameter vector
        data: Pandas data frame
    Returns:
        e: residual vector
    '''
    # determine number of entries in data frame
    n = len(data)

    # initialize matrix of residuals
    # rows: each row of Pandas data frame
    # columns: species CA, CB, CC
    e = np.zeros(n)

    # loop over experiments
    for i in data.exp.unique():

        # select the rows that correspond to the specific experiment number
        j = (data.exp == i)

        # determine experiment conditions
        CA0_ = float(data.CA0[j].mode())
        T_ = float(data.temp[j].mode())

        # determine experiment time
        # for newer versions of Pandas, use .to_numpy() instead
        t = data.time[j].to_numpy()

        CA, CB, CC = full_physics_model_2(theta,t,CA0_,T_)

        # only use observations of CB
        e[j] = CB - data.CB[j]

    return e

#-------------------------------------------------------------------
#KJ checkpoint in progress: dataframe is generated in reaction_genExpt_v3.ipynb

#data_temp = pd.read_csv()
#data = pd.DataFrame(['1','']
#e = regression_func_full_physics_model_2(theta, data)
#print('KJ checkpoint: nonlinear parameter estimation w/ full physics model')
#-------------------------------------------------------------------

def simple_physics_model(theta,t,CA0,T):
    '''Evaluate model A analytically

    Arugments:
        t: time, [hour], scalar or Numpy array
        theta: fitted parameters: A1, A2, E1, E2
        CA0: initial concentration, [mol/L], scalar or Numpy array
        T: temperature, [K], scalar or Numpy array

    Returns:
        CA, CB, CC: Concentrations at times t, [mol/L], three scalars or numpy arrays
    '''

    # units: [1/hr]
    k1 = kinetics(theta[0], theta[2], T)

    # units: [1/hr]
    k2 = kinetics(theta[1], theta[3], T)

    # units: [mol / L]
    CA = CA0 * np.exp(-k1*t);
    CB = k1*CA0/(k2-k1) * (np.exp(-k1*t) - np.exp(-k2*t));
    CC = CA0 - CA - CB;

    return CA, CB, CC

#------------------------------------------
#KJ Checkpoint:
#C_simple = simple_physics_model(theta, t, CA0, T)
#print('simple_physics_model working. C_simple=',C_simple)

#------------------------------------------

# nonlinear parameter estimation with full physics model
def regression_func_simple_physics_model(theta, data):
    '''
    Function to define regression function for least-squares fitting
    Arguments:
        theta: parameter vector
        data: Pandas data frame
    Returns:
        e: residual vector
    '''
    # determine number of entries in data frame
    n = len(data)

    # initialize matrix of residuals
    # rows: each row of Pandas data frame
    # columns: species CA, CB, CC
    e = np.zeros(n)

    # loop over experiments
    for i in data.exp.unique():

        # select the rows that correspond to the specific experiment number
        j = (data.exp == i)

        # determine experiment conditions
        CA0_ = float(data.CA0[j].mode())
        T_ = float(data.temp[j].mode())

        # determine experiment time
        # for newer versions of Pandas, use .to_numpy() instead
        t = data.time[j].to_numpy()

        CA, CB, CC = simple_physics_model(theta,t,CA0_,T_)

        # Only use CB measurements
        e[j] = CB - data.CB[j]

#    return e
#---------------------------------------------
#KJ Checkpoint
#figure out how to understand data which is a dataframe
#e = regression_func_simple_physics_model(theta, data)
#data = C_simple
#df = pd.DataFrame(data)
#print('line 223',df)


#---------------------------------------------

# TODO: Adapt this function to work with a Pandas dataframe instead or read in data from a text file.
# The GP model will not have a simple Python function to calculate the model predictions.
def plot_experiments_model(theta, data, model_func, fname=None, return_plot=False):
    ''' Plot the experimental and model predictions together...
    one plot per batch experiment

    Args:
        theta: model parameters
        data: Pandas data frame
        model_func: Python function to evaluate the model. Either 'full_physics_model` or `simple_physics_model`
        fname: string for prefix to filename to save the plot. if None, then do not save

    Returns:
        Nothing

    '''

    w_in = 3. # width in inches
    h_in = 4.0 # height in inches
    dpi_fig = 1200 # pixel density in dots per inch

    # loop over experiments
    for i in data.exp.unique():

        # delcare figure object
        fig, ax = plt.subplots(figsize=(w_in,h_in))
        # select the rows that correspond to the specific experiment number
        j = (data.exp == i)

        # determine experiment conditions
        CA0_ = float(data.CA0[j].mode())
        T_ = float(data.temp[j].mode())



        # Evaluate model
        ### BEGIN SOLUTION
        t = np.linspace(min(data.time[j]),max(data.time[j]),51)
        CA, CB, CC = model_func(theta,t,CA0_,T_)
        ### END SOLUTION

        # Plot model
        plt.plot(t, CA, label="$C_{A}$",linestyle="-",color="blue",linewidth=1)
        plt.plot(t, CB, label="$C_{B}$",linestyle="-.",color="green",linewidth=1)
        plt.plot(t, CC, label="$C_{C}$",linestyle="--",color="red",linewidth=1)

        # Plot data
        plt.plot(data.time[j], data.CA[j], marker='o',linestyle="",color="blue",label=str(),markersize=4)
        ### BEGIN SOLUTION
        plt.plot(data.time[j], data.CB[j], marker='s',linestyle="",color="green",label=str(),markersize=4)
        plt.plot(data.time[j], data.CC[j], marker='^',linestyle="",color="red",label=str(),markersize=4)
        ### END SOLUTION

        plt.xticks(np.arange(0.0,1.1,0.2))
        plt.xlim(0,1)
        # plt.yticks(np.arange(1.0,5.1,0.2))
        # plt.ylim(1,5)

        # Add "extras" to the plot
        plt.xlabel("Time (h)",fontsize=16,fontweight='bold')
        plt.ylabel(r"Conc. (mol l$\mathbf{^{-1}}$]",fontsize=16,fontweight='bold')

        plt_title = "Expt %s T=%.1f K,  CA0=%.1f mol l$^{-1}$"%(i+1,T_,CA0_)
        # plt.title(plt_title)
        # plt.title("Experiment "+str(i)+":  T="+str(T_)+" K,  CA0="+str(CA0_)+" mol/L")
        plt.legend()
        plt.grid()
        set_aspect_ratio(ax,1.0)

        if fname is not None:
            plt.savefig(fname+'{0}.png'.format(i),bbox_inches='tight',dpi=dpi_fig)
            plt.close()
        else:
            plt.show()

# TODO: Adapt this function to work with a Pandas dataframe instead or read in data from a text file.
# The GP model will not have a simple Python function to calculate the model predictions.
def verify_model(data1, data2, fname=None):
    ''' Overlay the experimental observations from two
    different datasets in one plot per experiment to verify
    that the datasets are identical. Used to verify the model
    parameters used for data generation

    Args:
        theta: model parameters
        data: Pandas data frame
        model_func: Python function to evaluate the model. Either 'full_physics_model` or `simple_physics_model`
        fname: string for prefix to filename to save the plot. if None, then do not save

    Returns:
        Nothing

    '''

    # loop over experiments
    for i in data1.exp.unique():

        # delcare figure object
        fig, ax = plt.subplots(figsize=(12,10))
        # select the rows that correspond to the specific experiment number
        j = (data1.exp == i)

        # determine experiment conditions
        CA0_ = float(data1.CA0[j].mode())
        T_ = float(data1.temp[j].mode())

        # Plot dataset 1
        plt.plot(data1.time[j], data1.CA[j], marker='o',markersize=16,linestyle="",color="blue",label="$C_{A} Set 1$")
        plt.plot(data1.time[j], data1.CB[j], marker='s',markersize=16,linestyle="",color="green",label="$C_{B} Set 1$")
        plt.plot(data1.time[j], data1.CC[j], marker='^',markersize=16,linestyle="",color="red",label="$C_{C} Set 1$")

        # Plot dataset 2
        plt.plot(data2.time[j], data2.CA[j], marker='o',markersize=10,linestyle="",color="fuchsia",label="$C_{A} Set 2$")
        plt.plot(data2.time[j], data2.CB[j], marker='s',markersize=10,linestyle="",color="fuchsia",label="$C_{B} Set 2$")
        plt.plot(data2.time[j], data2.CC[j], marker='^',markersize=10,linestyle="",color="fuchsia",label="$C_{C} Set 2$")

        # Add "extras" to the plot
        plt.xlabel("Time [hours]")
        plt.ylabel("Concentration [mol/L]")
        plt.title("Experiment "+str(i)+":  T="+str(T_)+" K,  CA0="+str(CA0_)+" mol/L")
        plt.legend()
        plt.grid()
        if fname is not None:
            plt.savefig(fname+'{0}.png'.format(i),bbox_inches='tight')
            plt.close()
        else:
            plt.show()
    




