##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 8
#
# Purpose  : Using the scipy.signal module to perform analysis of systems with rational polynomial transfer functions,coupled system of differential
#            equations, and also a linear electrical circuit which behaves like a low-pass filter.
# Author   : Neham Jain (EE19B084)
# Input    : No command line input is required
# Output   : The graphs are saved in a directory called plots

##**********************************************************************************************************************************************************************

#Libraries used in our program
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import scipy.signal as sp
import os
import warnings
warnings.filterwarnings("ignore")
from sympy import *

def sympy_to_lti(xpr, s=symbols('s')):
    """ Convert Sympy transfer function polynomial to Scipy LTI """

    num, den = simplify(xpr).as_numer_denom()  # returns the expressions
    p_num_den = poly(num, s), poly(den, s)
    num,den = p_num_den[0].all_coeffs(), p_num_den[1].all_coeffs()
    
    num=list(map(float, num))
    den=list(map(float, den))
    
    return sp.lti(num, den)

def mixed_freq_sinusoid(t):
    ''' Function that gives the sum of sinusoid with different frequencies '''
    return (np.sin(2000*np.pi*t)+np.cos(2000000*np.pi*t))*np.heaviside(t,0.5)

def damped_sinusoid(t,decay,freq):
    ''' Function that gives  a damped sinusoid '''

    return np.cos(freq*t)*np.exp(-decay*t) * (t>0)

def lowpass(R1=10e3,R2=10e3,C1=1e-9,C2=1e-9,G=1.586,Vi=1,s=symbols('s')):
    """Solve the given lowpass filter circuit for a given input Vi."""

    A=Matrix([[0,0,1,-1/G],
              [-1/(1+s*R2*C2),1,0,0],
              [0,-G,G,1],
              [-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    
    b=Matrix([0,0,0,-Vi/R1])    
    V=A.inv()*b
    
    return V

def highpass(R1=10e3,R3=10e3,C1=1e-9,C2=1e-9,G=1.586,Vi=1,s=symbols('s')):
    """Solve the given highpass filter circuit for a given input Vi."""
    
    A=Matrix([[0,-1,0,1/G],
        [s*C2*R3/(s*C2*R3+1),0,-1,0],
        [0,G,-G,1],
        [-1*s*C2-1/R1-s*C1,0,s*C2,1/R1]])

    b=Matrix([0,0,0,-Vi*s*C1])

    V=A.inv()*b
    
    return V

def get_bode_plot(H,title):
    ''' Plots the magnitude and phase response of a given transfer function '''

    w, mag, phase = H.bode()
    #Plotting Magnitude Response of Transfer Function 
    p = General_Plotter("Frequency in rad/s (log)","Magnitude in dB","Magnitude plot"+title)
    p.semilogx(w,mag)
    
    #Phase Response of Transfer Function
    p = General_Plotter("Frequency in rad/s (log)","Phase in degrees","Phase plot"+title)
    p.semilogx(w,phase)
    



class General_Plotter():
    ''' Class used for plotting different plots. Shortens the code by quite a bit'''
    
    fig_num=0   #Defined static variable for the figure number
    def __init__(self,xlabel,ylabel,title):
        ''' xlabel,ylabel,title are used in every graph''' 

        self.xlabel=xlabel
        self.ylabel=ylabel
        self.title=title
        self.fig=plt.figure(self.__class__.fig_num)
        self.__class__.fig_num+=1

    def general_funcs(self,ax):
        ''' General functions for every graph'''
        
        ax.set_ylabel(self.ylabel)
        ax.set_xlabel(self.xlabel)
        ax.set_title(self.title)
        #self.fig.show()
        self.fig.savefig("plots/"+self.title+".png")

    def general_plot(self,X,Y):
        ''' Line Plot'''
        axes=self.fig.add_subplot(111)
        axes.plot(X,Y)
        self.general_funcs(axes)
        
    
    def semilogx(self,X,Y):
        ''' Semilogx Plot ''' 
        
        axes=self.fig.add_subplot(111)
        axes.semilogx(X,Y)
        axes.grid()
        self.general_funcs(axes)

def main():
    
    #Getting Transfer Function for Lowpass circuit
    V_l = lowpass()[3]
    H1 = sympy_to_lti(V_l)
    
    #Plotting the Bode Plot of Transfer Function
    get_bode_plot(H1,title="1")
    
    #Step Response of the Lowpass Circuit
    t = np.linspace(0,0.001,1000)
    t,V_l_step = sp.step(H1,T=t)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Step Response of Low Filter Circuit")
    p.general_plot(t,V_l_step)

    #Step Response of the Highpass Circuit
    V_h = highpass()[3]
    H2 = sympy_to_lti(V_h)
    get_bode_plot(H2,title="2")
    t,V_h_step = sp.step(H2,T=t)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Step Response of High Filter Circuit")
    p.general_plot(t,V_h_step)
    
    #Plotting sum of sinusoids and the responses of circuits to the sum Of sinusoids
    t_ls = np.linspace(0,0.001,100000)
    t_hs = np.linspace(0,0.00001,100000)

    Vi = mixed_freq_sinusoid(t_ls)
    Vi_h = mixed_freq_sinusoid(t_hs)

    p=General_Plotter(r't$\rightarrow$',r'$V_{i}\rightarrow$',"Sum Of Sinusoids")
    p.general_plot(t_ls,Vi)
    
    t,Vo,_ = sp.lsim(H1,Vi,t_ls)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Response of Low Pass Filter Circuit to Sum of Sinusoids")
    p.general_plot(t,Vo)
    
    t,Vo,_ = sp.lsim(H2,Vi_h,t_hs)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Response of High Pass Filter Circuit to Sum of Sinusoids")
    p.general_plot(t,Vo)
    
    #Plotting the High Frequency damped Sinusoid and the response of circuits to this
    t_h = np.linspace(0,1e-3,100000)
    t_l = np.linspace(0,0.25,100000)
    
    Vi_h = damped_sinusoid(t_h,decay=3e3,freq=1e7)

    p=General_Plotter(r't$\rightarrow$',r'$V_{i}\rightarrow$',"High Frequency Decaying Sinusoid")
    p.general_plot(t_h,Vi_h)
    
    t,Vo,_ = sp.lsim(H1,Vi_h,t_h)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Response of Low Pass Filter Circuit to High Frequency Damped")
    p.general_plot(t,Vo)
    
    t,Vo,_ = sp.lsim(H2,Vi_h,t_h)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Response of High Pass Filter Circuit to High Frequency Damped")
    p.general_plot(t,Vo)
    
    #Plotting the Low Frequency damped Sinusoid and the response of circuits to this
    Vi_l = damped_sinusoid(t_l,decay=1e1,freq=1e3)

    p=General_Plotter(r't$\rightarrow$',r'$V_{i}\rightarrow$',"Low Frequency Decaying Sinusoid")
    p.general_plot(t_l,Vi_l)
    
    t,Vo,_ = sp.lsim(H1,Vi_l,t_l)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Response of Low Filter Circuit to Low Frequency Damped")
    p.general_plot(t,Vo)
    
    t,Vo,_ = sp.lsim(H2,Vi_l,t_l)
    p=General_Plotter(r't$\rightarrow$',r'$V_{o}\rightarrow$',"Response of High Filter Circuit to Low Frequency Damped")
    p.general_plot(t,Vo)
    
    print("The plots are saved at: ",os.getcwd()+"/plots/")


#if file is run directly
if __name__=="__main__":
    
    #Creating directory for storing plots
    os.makedirs("plots",exist_ok=True)
    #Running the main function
    print("Getting Results...")
    main()