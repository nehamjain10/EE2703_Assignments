##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 4
#
# Purpose  : To calculate the Fourier Series Coefficients for the given signals and verify the results using graphs.
# Author   : Neham Jain (EE19B084)
# Input    : Program requires no input 
# Output   : Plots different graphs as specified in the report
#
# NOTE: While plotting, we have tried to reduce the code redundancy as much as possible to minimize the size of the code
##**********************************************************************************************************************************************************************


#Importing the necessary libraries to be used in our program
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import os
import math

#Defining the exponential function
def exponential(x):
    return np.exp(x)    

#Defining the cos(cos(x)) function
def cos_cos(x):
    return np.cos(np.cos(x))


#extending exponential function to be periodic with given time period
def exp_extension(x):
    time_period=2*np.pi  
    return exponential(x%time_period)

#Defining the function f(x)*cos(kx) to calculate fourier coefficients
def mul_by_cos(x,func,k):
    return func(x)*np.cos(k*x)

#Defining the function f(x)*sin(kx) to calculate fourier coefficients
def mul_by_sin(x,func,k):
    return func(x)*np.sin(k*x)


#This function returns the first n Fourier Series Coefficients of the function f using the integration method 
def calculate_fourier_series_coefffs(n,function):
    a = np.zeros(n)         #List which stores all the coefficients
    a[0] = scipy.integrate.quad(function,0,2*math.pi)[0]/(2*math.pi) #Calculating a0
    for i in range(1,n):
        if(i%2==1):
            a[i] = scipy.integrate.quad(mul_by_cos,0,2*math.pi,args=(function,int(i/2)+1))[0]/math.pi  #Calculating an
        else:
            a[i] = scipy.integrate.quad(mul_by_sin,0,2*math.pi,args=(function,int(i/2)+1))[0]/math.pi  #Calculating bn
    return a 

#Plots the Fourier Coefficients in semilogy and loglog scale
def plot_fourier_coeffs(coeffs,type,ylabel,title,xlabel="Fourier Series Coefficients"):
    plt.figure()
    if type=="semilogy":
        plt.semilogy(np.abs(coeffs),'ro',label=r"Fourier Series Coefficients")
    if type=="loglog":
        plt.loglog(np.abs(coeffs),'ro',label=r"Fourier Series Coefficients")
    plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.savefig(f"plots/{title}.jpg")

#using lstsq for finding value of matrix c that best fits the equation Ac=b
def matrix_method(A,f,x):                         
    return scipy.linalg.lstsq(A,f(x))[0]            
    
#Function for plotting general functions
def general_func_plot(x_vals,y_vals,title,fmt,type="semilogy",xlabel="x",ylabel="y"):
    plt.figure()
    plt.grid(True)
    if type =="semilogy":
        plt.semilogy(x_vals,y_vals,fmt)
    elif type =='log':
        plt.loglog(x_vals,y_vals,fmt)
    elif type =="normal_scale":
        plt.plot(x_vals,y_vals,fmt)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(f"plots/{title}.jpg")

#Plotting aperiodic e^x and periodically extended e^x in log scale from x = -2pi to +4pi by splitting the x range into 10000 points from -2pi to 4pi 
def plot_exponential(x_vals):
    plt.figure(1) #figure-1
    plt.semilogy(x_vals,exponential(x_vals),'m',label='Aperiodic Exponential')
    plt.semilogy(x_vals,exp_extension(x_vals),'g',label='Periodically Extended Exponential')
    plt.grid(True)
    plt.ylabel(r'$log(e^{x})$',fontsize=15)
    plt.xlabel(r'$x$',fontsize=15)
    plt.legend()
    plt.savefig("plots/exponential.jpg")

#Plotting the coefficents obtained by matrix method and integration method and see deviation between the two approaches
def comparing_coeffs(coeffs_int,coeffs_mat,type,ylabel,title,xlabel="Magnitude Value"):
    plt.figure()
    if type=="semilogy":
        plt.semilogy(np.abs(coeffs_int),'go',label=r'Integration Approach')
        plt.semilogy(np.abs(coeffs_mat),'bo',label=r"Least Squares Approach")
    if type=="loglog":
        plt.loglog(np.abs(coeffs_int),'go',label=r'Integration Approach')
        plt.loglog(np.abs(coeffs_mat),'bo',label=r'Least Squares Approach')
    plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.savefig(f"plots/{title}.jpg")

#Function to plot the values of coefficients found by the 2 methods to see the deviations between the two functions
def plotting_convergence(fourier_func,f,x_vals,title,xlabel,ylabel):                 
    plt.figure()
    plt.semilogy(fourier_func, 'm', label = 'Fourier representation')
    plt.semilogy(f(x_vals), 'c', label = 'Original function')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.savefig(f"plots/{title}.jpg")    


#Main function from which all other functions are called
def main():

    #Create directory for storing plots
    os.makedirs("plots",exist_ok=True)
    
    #Question 1
    x_vals=np.linspace(-2*np.pi,4*np.pi,10000)
    plot_exponential(x_vals)        #Plotting exponential plot
    y_vals=cos_cos(x_vals)
    general_func_plot(x_vals,y_vals,"Cos_Cos_Plot","r-",type="normal_scale") #Plottting cos_cos plot

    #Question 2, Calculating Fourier Coefficients
    cos_cos_fsc = calculate_fourier_series_coefffs(51,cos_cos)
    exp_fsc = calculate_fourier_series_coefffs(51,exp_extension)
    
    #Question3, Plotting Fourier Series Coeffs in different scales
    plot_fourier_coeffs(cos_cos_fsc,type="semilogy",ylabel=r'$n$',title=r'Magnitude spectrum of coefficients in log scale for $cos(cos(x))$')
    plot_fourier_coeffs(cos_cos_fsc,type="loglog",ylabel=r'$n$',title=r'Magnitude spectrum of coefficients in loglog scale for $cos(cos(x))$')
    plot_fourier_coeffs(exp_fsc,type="loglog",ylabel=r'$n$',title=r'Magnitude spectrum of coefficients in loglog scale for $e^x$')
    plot_fourier_coeffs(exp_fsc,type="semilogy",ylabel=r'$n$',title=r'Magnitude spectrum of coefficients in log scale for $e^x$')  

    #Question4, Least Squares Approach
    x = np.linspace(0,2*np.pi,401)
    x=x[:-1]  
    A = np.zeros((400,51)) 
    A[:,0]=1 
    for k in range(1,26):
        A[:,2*k-1] = np.cos(k*x) 
        A[:,2*k] = np.sin(k*x) 
    
    #Question5, Finding the best fit of coefficients using matrix method
    coeffs_exp_lstsq = matrix_method(A,exponential,x)       
    coeffs_cos_cos_lstsq = matrix_method(A,cos_cos,x) 


    #Question6, Comparing deviation of coefficients between the two methods
    comparing_coeffs(exp_fsc,coeffs_exp_lstsq,type="semilogy",ylabel=r'$n$',title=r'Comparing magnitude spectrum of coefficients in log scale for $e^x$',xlabel="Magnitude Value")
    comparing_coeffs(exp_fsc,coeffs_exp_lstsq,type="loglog",ylabel=r'$n$',title=r'Comparing magnitude spectrum of coefficients in loglog scale for $e^x$',xlabel="Magnitude Value")
    comparing_coeffs(cos_cos_fsc,coeffs_cos_cos_lstsq,type="semilogy",ylabel=r'$n$',title=r'Comparing magnitude spectrum of coefficients in log scale for $cos(cos(x))$',xlabel="Magnitude Value")
    comparing_coeffs(cos_cos_fsc,coeffs_cos_cos_lstsq,type="loglog",ylabel=r'$n$',title=r'Comparing magnitude spectrum of coefficients in loglog scale for $cos(cos(x))$',xlabel="Magnitude Value")

    max_error_f1 = np.max(np.abs(exp_fsc - coeffs_exp_lstsq))
    max_error_f2 = np.max(np.abs(cos_cos_fsc - coeffs_cos_cos_lstsq))

    print("Maximum error for coefficients of e^x is ",max_error_f1)
    print("Maximum error for coefficients of cos(cos(x)) is ",max_error_f2)
    
    #Question7, Comparing deviation of functions between the two methods 
    fourier_func_exp = np.matmul(A,coeffs_exp_lstsq)
    plotting_convergence(fourier_func_exp,exp_extension,x,title=r'Convergence of Fourier Series representation to actual function for e^x',xlabel=r'$x$',ylabel=r'Value in log scale')    
    fourier_func_cos = np.matmul(A,coeffs_cos_cos_lstsq)
    plotting_convergence(fourier_func_cos,cos_cos,x,title=r'Convergence of Fourier Series representation to actual function for cos(cos(x))',xlabel=r'$x$',ylabel=r'Value in log scale')

main()