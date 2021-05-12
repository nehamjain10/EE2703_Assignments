##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 7
#
# Purpose  : To analyse filter circuits using the Laplace Transform of the impulse response, utilising the Symbolic Algebra capabilities of Python.
# Author   : Neham Jain (EE19B084)
# Input    : No command line input is required
# Output   : The graphs are saved in a new directory called plots

##**********************************************************************************************************************************************************************

#Libraries used in our program
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as sp
import os
import sympy

def input_signal_laplace(freq,decay):
    """Laplace Transform of the Input Signal given in Question 1"""
    
    n = np.poly1d([1,decay])
    d = n*n+freq**2
    return n,d

def general_transfer_fcn(wn=1.5,zeta=0,gain=1/2.25):
    """General transfer function for a second order system"""
    
    n = np.poly1d([wn**2*gain])
    d = np.poly1d([1,2*wn*zeta,wn**2])
    return n,d

def lti_solver(decay,freq=1.5):
    """Find the response to the given system to a decaying cosine."""
    
    input_numerator, input_denominator = input_signal_laplace(freq,decay=decay)
    transfer_numerator, transfer_denominator = general_transfer_fcn()

    output_numerator,output_denominator = input_numerator*transfer_numerator, input_denominator*transfer_denominator
    out_s = sp.lti(output_numerator.coeffs, output_denominator.coeffs)

    t = np.linspace(0,50,1000)
    
    return sp.impulse(out_s,None,t)

def input_signal_time(t,decay=0.5,freq=1.5):
    """Input Signal in Time Domain of the signal given in Question 1"""

    u_t = 1*(t>0)
    return np.cos(freq*t)*np.exp(-decay*t) * u_t

def cosines(t,w1=1e3,w2=1e6):
    """Two cosines of different frequencies"""
    u_t = 1*(t>0)
    return (np.cos(w1*t)-np.cos(w2*t)) * u_t


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

    def general_plot(self,X,Y,legend_txt=None):
        ''' Line Plot'''
        axes=self.fig.add_subplot(111)
        axes.plot(X,Y)
        if legend_txt is not None:
            axes.legend(labels=legend_txt)
        self.general_funcs(axes)
        
    
    def semilogx(self,X,Y):
        ''' Semilogx Plot ''' 
        axes=self.fig.add_subplot(111)
        axes.semilogx(X,Y)
        axes.grid()
        self.general_funcs(axes)



def main():
    
    #Question 1
    X1,Y1=lti_solver(0.5)
    p1=General_Plotter(r"$t$",r"$x$",r"System response with decay of 0.5")
    p1.general_plot(X1,Y1)

    #Question 2
    X2,Y2=lti_solver(0.05)
    p2=General_Plotter(r"$t$",r"$x$",r"System response with decay of 0.05")
    p2.general_plot(X2,Y2)

    #Question 3
    transfer_fcn=general_transfer_fcn(wn=1.5,zeta=0,gain=1/2.25)
    outs = []

    #List of Frequencies to iterate over
    freqs = np.linspace(1.4,1.6,5)
    t = np.linspace(0,70,1000)
    for freq in freqs:
        #Solving
        t,y,_ = sp.lsim(transfer_fcn,input_signal_time(t,decay=0.05,freq=freq),t)
        #Storing in list outs
        outs.append(y)
    
    #Plotting System Response with variation of input frequency    
    p3 = General_Plotter(r"$t$",r"$x$",r"System responses with variation of input frequency")
    p3.general_plot(t,np.array(outs).T,["Freq = ${:.2f}$".format(f) for f in freqs])

    #Bode Plot
    w,S,phi=sp.lti(*transfer_fcn).bode()

    #Plotting Magnitude Response of Transfer Function 
    p4 = General_Plotter("Frequency in rad/s (log)","Magnitude in dB","Magnitude plot")
    p4.semilogx(w,S)
    
    #Phase Response of Transfer Function
    p5 = General_Plotter("Frequency in rad/s (log)","Phase in degrees","Phase plot")
    p5.semilogx(w,phi)
    
    #Question 4
    X_s = sp.lti([1,0,2],[1,0,3,0])
    Y_s = sp.lti([2],[1,0,3,0])

    #Plotting X and Y impulse
    t = np.linspace(0,20,1000)
    t, x = sp.impulse(X_s,None,t)
    t, y = sp.impulse(Y_s,None,t)

    #Plotting Displacement of Coupled System
    p6 = General_Plotter(r"$t$","Displacement","Responses of coupled system")
    p6.general_plot(t,np.array([x,y]).T,legend_txt=[r"$x(t)$", r"$y(t)$"])

    #Question 5
    #Find the transfer function of the given circuit
    R = 100
    L = 1e-6
    C = 1e-6

    wn = 1/np.sqrt(L*C) #Natural Frequency
    Q = 1/R * np.sqrt(L/C)  #Quality Factor
    zeta = 1/(2*Q)  #Damping Constant

    #Transfer function
    n,d = general_transfer_fcn(gain=1,wn=wn,zeta=zeta)

    #Make System
    H = sp.lti(n,d)

    #Get Bode Plots
    w,S,phi=H.bode()

    #Plotting Magnitude Response of Transfer Function 
    p7 = General_Plotter("Frequency in rad/s (log)","Magnitude in dB","Magnitude plot 2")
    p7.semilogx(w,S)
    
    #Phase Response of Transfer Function
    p8 = General_Plotter("Frequency in rad/s (log)","Phase in degrees","Phase plot 2")
    p8.semilogx(w,phi)

    #Question 8
    #Transient Response
    t1=np.linspace(0,30e-6,1000)
    t1,y1,_ = sp.lsim(H,cosines(t1),t1)

    #Steady State Response
    t2=np.linspace(0,10e-3,1000)
    t2,y2,_ = sp.lsim(H,cosines(t2),t2)

    #Plotting Transient Response
    p9=General_Plotter(r"$t$ (sec)",r"$v_0(t)$",r"Response for 30 micro seconds")
    p9.general_plot(t1,y1)

    #Plotting Steady State Response 
    p10=General_Plotter(r"$t$ (sec)",r"$v_0(t)$",r"Response for 10 msec")
    p10.general_plot(t2,y2)
    
    print("The plots are saved at: ",os.getcwd()+"/plots/")


#if file is run directly
if __name__=="__main__":
    
    #Creating directory for storing plots
    os.makedirs("plots",exist_ok=True)
    #Running the main function
    main()