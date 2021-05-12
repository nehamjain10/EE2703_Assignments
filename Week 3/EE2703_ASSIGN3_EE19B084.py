##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 3

# Purpose  : To determine how the amount of noise affects the quality of estimation of the function
# Author   : Neham Jain (EE19B084)
# Input    : Program takes in the columns of fitting.dat as input
# Output   : Plots different graphs as specified in the report

##**********************************************************************************************************************************************************************

#Importing the required libraries to be used (we have not used pylab as it's use is discouraged due to NameSpace polluting)

from pylab import *
import scipy.special as sp
import warnings
import sys
import os
warnings.filterwarnings("ignore")

#Function for calculating g(t;A,B)
def g(t,A,B):
    return(A*sp.jn(2,t)+B*t) # sp.jn denotes the bessel function

#Function for calculating MSE
def calculate_mse(predictions, targets):
    return ((predictions - targets) ** 2).mean()

#Function for constructing the M matrix
def construct_M(time):
    return c_[sp.jn(2,time),time]

#Estimating the value of A and B
def find_AB(M,b):
    return linalg.lstsq(M,b) 

#Helper Function for Question 3 and 4
def plot_fitting_cols(time,vals,ground_truths,scl):
    vals=c_[vals,ground_truths]
    figure(0)
    plot(time,vals)
    xlabel(r"Time $\rightarrow$")
    ylabel(r"$f(t)+Noise \rightarrow$")
    title("Plotting fitting.dat and ground truth values")
    scl=[f"Standard Deviation={round(i,4)}" for i in scl]
    legend(list(scl)+["True Value"])
    grid(True)
    #show()
    savefig("plots/Ques3_4.jpg")
    close()

#Helper function for plotting the error bars
def plot_errorbars(time,values,ground_truth):
    figure(1)
    std_dev = std(values[:,0]-ground_truth)
    xlabel(r"Time $\rightarrow$")
    ylabel(r"f(t) $\rightarrow$")
    title(f"Data Points for stddev={round(std_dev,3)} along with exact function")
    errorbar(time[::5],values[:,0][::5],std_dev,fmt='ro')
    plot(time,ground_truth)
    legend(["f(t)","Error Bar"])
    #show()
    savefig("plots/Ques5.jpg")
    close()

#Helper function for plotting the contour plot
def plot_contours(values,time):
    figure(2)
    A_range=arange(0,2.1,0.1)
    B_range=arange(-0.2,0.01,0.01)
    epsilon_matrix = zeros((len(A_range),len(B_range)))
    for count_A,A in enumerate(A_range):
        for count_B,B in enumerate(B_range):
            epsilon_matrix[count_A][count_B] = calculate_mse(values[:,0],g(time,A,B))
    #print(epsilon_matrix)
    contour_obj = contour(A_range,B_range,epsilon_matrix,levels=arange(0,20*0.025,0.025),cmap="magma")
    clabel(contour_obj,contour_obj.levels[0:5],inline=True,fontsize=10)
    title('Contour Plot of Error')
    xlabel(r'A $\rightarrow$',size=12)
    ylabel(r'B $\rightarrow$',size=12)
    plot(1.05, -0.105,'ro', label = 'Exact Value')
    annotate("Exact Value",xy = [0.8,-0.100])
    savefig("plots/Ques8.jpg")
    #show()
    close()

#Helper function for plotting the contour plot
def plot_variation_of_error(scl,error_a,error_b):
    figure(3)
    plot(scl,error_a,'r--')
    scatter(scl,error_a)
    plot(scl,error_b, 'b--')
    scatter(scl,error_b)
    legend(["Aerr","Berr"])
    title("Variation Of Error with Noise")
    xlabel(r'$\sigma_{n}\rightarrow$',size=10)
    ylabel(r'Absolute Error $\rightarrow$',size=10)
    savefig("plots/Ques10.jpg")
    close()

#Helper function for plotting the log variation of error and the stem plot
def plot_log_variation_of_error(scl,error_a,error_b):
    figure(5)
    stem(scl,error_a,'ro')
    loglog(scl,error_a,'ro')
    loglog(scl,error_b,'bo')
    stem(scl,(error_b),'bo')
    xlabel(r'$\sigma_{n}\rightarrow$',fontsize=15)
    ylabel(r'Error$\rightarrow$',fontsize=15)   
    legend(["Aerr","Berr"])
    title("Stem plot showing variation of error with standard deviation")
    savefig("plots/Ques11.jpg")
    #show()
    close()

def main():

    #Create directory for storing plots
    os.makedirs("plots",exist_ok=True)

    #Question 2 (loading the data)
    try:
        data=loadtxt("fitting.dat")
    except:
        print("fitting.dat not found. Please run generate_data.py again")
        sys.exit(0)

    #Parsing the data
    time_values=data[:,0]
    func_values=data[:,1:]

    #Calculating ground truth values for Question 6
    ground_truth = g(time_values,1.05,-0.105)

    #Plotting required graphs

    #Question 3 and 4
    scl=logspace(-1,-3,9) # noise stdev
    plot_fitting_cols(time_values,func_values,ground_truth,scl)

    #Question 5
    plot_errorbars(time_values,func_values,ground_truth)

    #Question 6
    M=construct_M(time_values)
    AB=np.array([1.05,-0.105])
    #Confirms whether M*AB and the function g(t,A,B) are equal or not
    print("The mean squared error between M*AB and g(t,A,B) is: ",calculate_mse(matmul(M,AB),g(time_values,1.05,-0.105)))

    #Question 7 and 8
    plot_contours(func_values,time_values)

    #Question 9
    AB_gt_values=find_AB(M,ground_truth)[0] #The ground truth values of A and B 
    print("The best estimate of A and B for ground truth values are: ", round(AB_gt_values[0],3),round(AB_gt_values[1],3))
    
    #Question 10
    a_error=zeros(9)
    b_error=zeros(9)
    total_error=zeros(9)
    
    #Iterating through all the columns of fitting.dat file and calculating the absolute error between true values of A and B vs the predicted value
    for i in range(9):
        temp = find_AB(M,func_values[:,i])
        prediction=temp[0]
        diff=temp[1]
        a_error[i],b_error[i] = abs(prediction[0]-AB_gt_values[0]),abs(prediction[1]-AB_gt_values[1])
        total_error[i] = diff
    
    plot_variation_of_error(scl,a_error,b_error)
    
    #Question 11
    plot_log_variation_of_error(scl,a_error,b_error)

    print("Plots are saved in following directory: ",os.getcwd()+"/plots/")

#Run the main function
main()