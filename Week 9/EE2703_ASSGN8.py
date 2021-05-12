##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 8
#
# Purpose  : Analysing the frequency spectrum of signals using numpy.fft module. and approximate the CTFT of a Gaussian
# Author   : Neham Jain (EE19B084)
# Input    : No command line input is required
# Output   : The graphs are saved in a directory called plots

##**********************************************************************************************************************************************************************

#Libraries used in our program
import numpy as np 
import matplotlib
from pylab import *
import matplotlib.pyplot as plt 
import os

#Setting Default Parameters for the Plots
rcParams['figure.figsize'] = 12,8
rcParams['axes.grid'] = True
rcParams['font.size'] = 14


class General_Plotter():
    ''' Class used for plotting different plots. Shortens the code by quite a bit'''
    
    fig_num=0   #Defined static variable for the figure number
    def __init__(self,xlabel1,ylabel1,ylabel2,title1,save_name):
        ''' xlabel,ylabel,title are used in every graph''' 

        self.xlabel1 = xlabel1
        self.ylabel1 = ylabel1
        
        self.xlabel2 = xlabel1
        self.ylabel2=ylabel2
        
        self.title1=title1

        self.save_name=save_name
        self.fig=plt.figure(self.__class__.fig_num)
        self.__class__.fig_num+=1

    def general_funcs(self,ax,xlim,ylim):
        ''' General functions for every graph'''
        
        ax[0].set_ylabel(self.ylabel1)
        ax[0].set_xlabel(self.xlabel1)
        ax[0].set_title(self.title1)
        ax[1].set_ylabel(self.ylabel2)
        ax[1].set_xlabel(self.xlabel2)
        if xlim is not None:
            ax[0].set_xlim(-xlim,xlim)
            ax[1].set_xlim(-xlim,xlim)
        if ylim is not None:
            ax[1].set_ylim(-ylim,ylim)

        plt.grid(True)
        plt.tight_layout()
        #plt.show()
        self.fig.savefig("plots/"+self.save_name+".png")

    def plot_fft(self,w,Y,scatter=False,ro=False,xlim=None,ylim=None):
        ''' Helper Function for plotting the fft of a given signal'''
        axes=self.fig.subplots(2,1)
        axes[0].plot(w,abs(Y),lw=2)

        if scatter==True:
            phi=angle(Y)
            phi[where(abs(Y)<1e-5)]=0
            if ro==True:
                axes[1].plot(w,phi,'ro')

            ii=where(abs(Y)>1e-5)[0]
            axes[1].plot(w[ii],phi[ii],'go')
        else:    
            axes[1].plot(angle(Y),lw=2)
        self.general_funcs(axes,xlim,ylim)


def estimate_dft(func_name,x_start,x_end,steps,xlim1,title1,title2,ylabel1,ylabel2,xlabel1,savename,ro=True,ylim=None):
    """
    Estimate the DFT of a given signal using the fft module of the numpy
    Phase of points of magnitude having very less is taken as zero
    """
    
    sampling_rate = steps/(x_end-x_start)
    
    x=linspace(x_start,x_end,steps+1)[:-1]
    y = func_name(x)
    
    Y=fftshift(fft(y))/float(steps)
    w=sampling_rate*(linspace(-pi,pi,steps+1)[:-1])
    
    p1=General_Plotter(xlabel1,ylabel1,ylabel2,title1,title2,savename)
    p1.plot_fft(w,Y,scatter=True,ro=ro,xlim=xlim1,ylim=ylim)

def estimateCTFT(func, tol=1e-6,time_samples=128, true_func=None,func_name=None, wlim=None, scatter_size=40):
    """
    Estimate the continuous time Fourier Transform of the given function
    by finding the DFT of a sampled window of the function. The magnitude and
    phase of the estimate are also plotted.
    """
    T = 8*pi
    N = time_samples
    Xold = 0
    error = tol+1
    iters=0
    
    while error>tol:
        
        delta_t = T/N # time resolution
        delta_w = 2*pi/T # frequency resolution

        W = N*delta_w # total frequency window size

        t = linspace(-T/2,T/2,N+1)[:-1] # time points
        w = linspace(-W/2,W/2,N+1)[:-1] # freq points

        x = func(t)

        # Find DFT and Normalize
        # note that ifftshift is used to prevent artifacts in the
        # phase of the result due to time domain shifting
        X = delta_t/(2*pi) * fftshift(fft(ifftshift(x)))
        
        error = sum(abs(X[::2]-Xold))
        
        Xold = X
        N *= 2 # number of samples
        T *= 2 # total time window size
        iters+=1
        
    print("\nDTFT Approximation Results:")    
    print("Estimated error after {} iterations: {}".format(iters, error))
    print("Time range : ({:.4f}, {:.4f})".format(-T/2,T/2))
    print("Time resolution : {:.4f}".format(delta_t))
    print("Frequency resolution : {:.4f}".format(delta_w))
        
    if true_func != None:
        true_error = sum(abs(X-true_func(w)))
        print("True error: {}".format(true_error))
    
    mag = abs(X)
    ph = angle(X)
    ph[where(mag<tol)]=0
    
    #Magnitude
    p1=General_Plotter("Frequency in rad/s","Magnitude",r"Phase",r"Magnitude of CFT Estimate of $e^{\frac{-t^2}{2}}$",r"Phase of CFT Estimate","CFT_Estimate")
    p1.plot_fft(w,mag*np.exp(1j*ph),scatter=True,ro=True,xlim=wlim)
    
    X_ = true_func(w)    
    mag = abs(X_)
    ph = angle(X_)
    ph[where(mag<tol)]=0
    
    p2=General_Plotter("Frequency in rad/s","Magnitude",r"Phase",r"Magnitude of True CFT of $e^{\frac{-t^2}{2}}$",r"Phase of True CFT","True_CFT")
    p2.plot_fft(w,mag*np.exp(1j*ph),scatter=True,ro=True,xlim=wlim)
    

#Definitions of the functions whose fft needs to be calculated

def sin5(t): 
    return sin(5*t)

def amp_mod(t):
    return (1+0.1*cos(t))*cos(10*t)

def sin3(t):
    return sin(t)**3

def cos3(t):
    return cos(t)**3

def freq_mod(t):
    return cos(20*t + 5*cos(t))

#defining gaussian and its expected CTFT 
def gauss(x):
    return exp(-0.5*x**2)

def expectedgauss(w):
    return 1/sqrt(2*pi) * exp(-w**2/2)

def main():
    
    #Assignment_Examples

    #Example 1
    x=rand(100)
    X=fft(x)
    y=ifft(X)
    c_[x,y]
    print ("The Absolute Maximum Error is ",abs(x-y).max())

    #Example 2
    x=linspace(0,2*pi,128)
    y=sin(5*x)
    Y=fft(y) #finding fft

    p1=General_Plotter("Frequency in rad/s","Magnitude",r"Phase",r"Magnitude of DFT of $sin(5x)$",r"Phase of DFT of $sin(5x)$","fft_unshifted_sin(5x)")
    p1.plot_fft(arange(len(Y)),Y)

    #Example 3
    estimate_dft(sin5,0,2*pi,128,10,r"Spectrum of $\sin(5t)$",r"Phase of $Y$","Magnitude","Phase","Frequency in rad/s","sin(5x)_shifted",ro = True)

    #Example 4
    estimate_dft(amp_mod,0,2*pi,128,15,r"Spectrum of $(1+0.1*cos(t))*cos(10t)$",r"Phase of $Y$","Magnitude","Phase","Frequency in rad/s","1+cos(0.1t)",ro = True,ylim=1.5)

    #Example 5
    estimate_dft(amp_mod,-4*pi,4*pi,512,15,r"Spectrum of $(1+0.1*cos(t))*cos(10t)$",r"Phase of $Y$","Magnitude","Phase","Frequency in rad/s","1+cos(0.1t)_stretched",ro = True,ylim=1.5)

    #Assignment Questions

    #DTFT of sin^3(t)
    estimate_dft(sin3,-4*pi,4*pi,512,15,r"Spectrum of $sin^3(t)$",r"Phase of $Y$","Magnitude","Phase","Frequency in rad/s","sin3",ro = True)

    #DTFT of cos^3(t)
    estimate_dft(cos3,-4*pi,4*pi,512,15,r"Spectrum of $cos^3(t)$",r"Phase of $Y$","Magnitude","Phase","Frequency in rad/s","cos3",ro = True,ylim=1.5)

    #DTFT of Frequency Modulated Signal
    estimate_dft(freq_mod,-4*pi,4*pi,512,40,r"Spectrum of $cos(20t + 5cos(t))$",r"Phase of $Y$","Magnitude","Phase","Frequency in rad/s","freq_mod",ro = False)

    #Estimate CTFT
    estimateCTFT(gauss,true_func=expectedgauss,wlim=10,func_name=r"$e^{\frac{-t^2}{2}}$")
    
    print("The plots are saved at: ",os.getcwd()+"/plots/")



#if file is run directly
if __name__=="__main__":
    
    #Creating directory for storing plots
    os.makedirs("plots",exist_ok=True)
    #Running the main function
    print("Getting Results...")
    main()
