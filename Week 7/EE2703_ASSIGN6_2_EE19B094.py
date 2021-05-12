'''***************************************************************
PURPOSE : EE2703 - Assignment 6-2
AUTHOR : Manvar Nisharg (EE19B094)
INPUT : NULL
OUTPUT : 12 plots (2--Q1,2)(5--Q3)(2--Q2)(1--Q5)(2--Q6)
***************************************************************'''


import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy

class General_Plotter():
    ''' Class used for plotting different plots. Shortens the code by quite a bit'''
    
    def __init__(self,xlabel,ylabel,title,legend,fig_num=0):
        ''' xlabel,ylabel,title are used in every graph''' 

        self.xlabel=xlabel
        self.ylabel=ylabel
        self.title=title
        self.legend=legend
        self.fig=plt.figure(fig_num)
    
    def general_funcs(self,ax):
        ''' General functions for every graph'''
        
        ax.set_ylabel(self.ylabel)
        ax.set_xlabel(self.xlabel)
        ax.set_title(self.title)
        #legend = ax.legend()
        plt.show()
        #self.fig.savefig(self.title+".png")


    def line_plot(self,X,Y):
        axes=self.fig.add_subplot(111)
        axes.plot(X,Y)
        self.general_funcs(axes)


#Q1
def Q_1(freq,decay_constant):
	#Function to get laplace equation for input for given values of frequency and decay constant
	def input_laplace(freq,decay_constant):
		num = np.poly1d([1,decay_constant])
		den = np.poly1d([1,2*decay_constant,decay_constant*decay_constant+freq*freq])
		return num,den

	#Input function
	input_num,input_den = input_laplace(freq,decay_constant)

	#Output function by multiplying transfer function and input function
	output_laplace_num = input_num
	output_laplace_den = np.polymul(input_den,[1,0,2.25])
	output_laplace = sp.lti(output_laplace_num,output_laplace_den)

	#Get time domain output for time = 0 to 5001
	time,output_time = sp.impulse(output_laplace,None,np.linspace(0,50,5001))

	#Plot graph
	g1 = General_Plotter("Time$\longrightarrow$","Position$\longrightarrow$","Q1: Decay = %s"%decay_constant,[])
	g1.line_plot(time,output_time)


Q_1(1.5,0.5)	#Decay constant = 0.5
Q_1(1.5,0.05)	#Decay constant = 0.05




#Q3
def Q_3():
	#Helper function which returns input function values at given time stamps
	def input_time(freq,decay_constant,time):
		cos = np.cos(freq*time)
		exp = np.multiply(np.exp(-decay_constant*time),np.heaviside(time,0.5))
		return np.multiply(cos,exp)

	#Looping through different values of frequency and plotting the output time domain simulation
	for freq in np.arange(1.4,1.6,0.05):
		transfer_function = sp.lti([1],[1,0,2.25])
		time_stamps = np.linspace(0,100,1000)

		#Simulating
		_,response,_ = sp.lsim(transfer_function,input_time(freq,0.05,time_stamps),time_stamps)
		g1 = General_Plotter("Time$\longrightarrow$","Position$\longrightarrow$","Q3: Frequency = %s"%freq,[])
		g1.line_plot(time_stamps,response)

Q_3()

#Q4
def Q_4():
	#Laplace function for x and y obtained by decopling the equations
	laplace_function_x = sp.lti(np.poly1d([1,0,2]),np.poly1d([1,0,3,0]))
	laplace_function_y = sp.lti(np.poly1d([2]),np.poly1d([1,0,3,0]))
	time_stamps = np.linspace(0,20,1000)

	#Obtaining the output in time domain for each laplace function
	response_x = sp.impulse(laplace_function_x,None,time_stamps)
	response_y = sp.impulse(laplace_function_y,None,time_stamps)

	g1 = General_Plotter("Time$\longrightarrow$","X$\longrightarrow$","Q4: X vs time",[])
	g1.line_plot(response_x[0],response_x[1])

	g1 = General_Plotter("Time$\longrightarrow$","Y$\longrightarrow$","Q4: Y vs time",[])
	g1.line_plot(response_y[0],response_y[1])

Q_4()

#Q5
def Q_5():
	#Transfer function for the system
	transfer_function = sp.lti(np.poly1d([1000000]),np.poly1d([0.000001,100,1000000]))

	#Obtaining and plotting the bode plot of the transfer function
	freq,magnitude,phase = transfer_function.bode()

	plt.subplot(2,1,1)
	plt.semilogx(freq,magnitude)
	plt.ylabel(r'$|H(s)|$')
	plt.subplot(2,1,2)
	plt.semilogx(freq,phase)
	plt.ylabel(r'$\angle(H(s))$')
	plt.suptitle("Q5: Bode plot of two port network")
	plt.show()
	return transfer_function

transfer_function = Q_5()

#Q6
def Q_6(transfer_function):
	#Simulating and plotting Q5 for given input for time uptill 30us 
	time = np.linspace(0,30*0.000001,1000)
	vi = np.multiply(np.cos(1000*time)-np.cos(1000000*time),np.heaviside(time,0.5))
	_,output_time,_ = sp.lsim(transfer_function,vi,time)
	g1 = General_Plotter("Time$\longrightarrow$","Voltage$\longrightarrow$","Q6: Voltage uptill 30$\mu$s",[])
	g1.line_plot(time,output_time)
	
	#Simulating and plotting Q5 for given input for time uptill 1ms
	time = np.linspace(0,10*0.001,100000)
	vi = np.multiply(np.cos(1000*time)-np.cos(1000000*time),np.heaviside(time,0.5))
	_,output_time,_ = sp.lsim(transfer_function,vi,time)
	g1 = General_Plotter("Time$\longrightarrow$","Voltage$\longrightarrow$","Q6: Voltage uptill 10ms",[])
	g1.line_plot(time,output_time)

Q_6(transfer_function)
