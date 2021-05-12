import sympy
import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy
import warnings
warnings.filterwarnings("ignore")
sympy.init_session


#Plotter class to plot graphs
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
		ax.grid(True)
		#self.fig.savefig(self.title+".png")

	def show(self,show_legend=0):
		if(show_legend):
			axes=self.fig.add_subplot(111)
			axes.legend()
		plt.show()


	def plot_loglog(self,X,Y,label_name=""):
		axes=self.fig.add_subplot(111)
		axes.loglog(X,Y,label=label_name)
		self.general_funcs(axes)

	def plot_points(self,X,Y,label_name=""):
		axes=self.fig.add_subplot(111)
		axes.plot(X,Y,label=label_name)
		self.general_funcs(axes)


#Convert transfer function in sympy to scipy
def convert_sympy_to_scipy(H):
	H = sympy.simplify(H)
	n,d = sympy.fraction(H)		#Get numerator and denominator from H
	num,den = sympy.Poly(n,s), sympy.Poly(d,s)	#Convert them into polynomial in 's'
	num,den = num.all_coeffs(), den.all_coeffs()	#Get coefficients of polynomials
	num,den = [float(f) for f in num], [float(f) for f in den]	#Store them in list
	return num,den


#Calculate step response to given transfer function in sympy
def step_response(H):
	num,den = convert_sympy_to_scipy(H)		#Sympy to scipy
	den.append(0)							#Multiply H by 1/s for step input
	H = sp.lti(num,den)
	t,x = sp.impulse(H,T = np.linspace(0,1e-3,100000))	#Calculate Impulse response
	return t,x


#Function to calculate output corresponding to given transfer function and input function
def general_input(H,input_time,max_time=1e-3):
	num,den = convert_sympy_to_scipy(H)		#Sympy to scipy
	H = sp.lti(num,den)
	t = np.linspace(0,max_time,100000)		#Time range
	t,y,svec = sp.lsim(H,input_time(t),t)	#Calculate output
	return t,y

#Sum of sinusoids of two frequencies
def sum_of_sinuosids(t):
	return (np.sin(2000*np.pi*t)+np.cos(2e6*np.pi*t))

#High frequency damped sinusoid
def damped1(t,decay=3e3,freq=1e7):
	return np.cos(freq*t)*np.exp(-decay*t) * (t>0)

#Low frequency damped sinusoid
def damped2(t,decay=1e1,freq=1e3):
	return np.cos(freq*t)*np.exp(-decay*t) * (t>0)


#KCL for low pass filter
def low_pass_filter(R1=10e3,R2=10e3,C1=1e-9,C2=1e-9,G=1.586,Vi=1):
	s = sympy.symbols('s')
	A = sympy.Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
	b = sympy.Matrix([0,0,0,-Vi/R1])
	V = A.inv()*b
	return (A,b,V)

#KCL for high pass filter
def high_pass_filter(R1=10e3,R3=10e3,C1=1e-9,C2=1e-9,G=1.586,Vi=1):
	s = sympy.symbols('s')
	A = sympy.Matrix([[0,-1,0,1/G], [s*C2*R3/(s*C2*R3+1),0,-1,0],[0,G,-G,1], [-s*C2-1/R1-s*C1,0,s*C2,1/R1]])
	b = sympy.Matrix([0,0,0,-Vi*s*C1])
	V = A.inv()*b
	return (A,b,V)

#Frequency sweep range
s = sympy.symbols('s')
w = np.logspace(0,8,801)
ss = 1j*w


#Bode Plots
#********************************************
A,b,V = low_pass_filter()	#Solve for low pass filter
H_lowpass = V[3]
hf = sympy.lambdify(s,H_lowpass,"numpy")	#Convert transfer function to python function
v = hf(ss)		#Do freq sweep and plot magnitude response
g1 = General_Plotter("frequency$\longrightarrow$","|H(s)|$\longrightarrow$","Bode plot of low pass filter",[],1)
g1.plot_loglog(w,abs(v))
g1.show()


A,b,V = high_pass_filter()	#Solve for high pass filter
H_highpass = V[3]
hf = sympy.lambdify(s,H_highpass,"numpy")		#Convert transfer function to python function
v = hf(ss)		#Do freq sweep and plot magnitude response
g1 = General_Plotter("frequency$\longrightarrow$","|H(s)|$\longrightarrow$","Bode plot of high pass filter",[],2)
g1.plot_loglog(w,abs(v))
g1.show()
#********************************************



#STEP RESPONSES
#********************************************
t,x = step_response(H_lowpass)	#Get stepresponse for low pass filter and plot
g1 = General_Plotter("time$\longrightarrow$","$V_o(t)\longrightarrow$","Step response for low pass filter",[],3)
g1.plot_points(t,x)
g1.show()

t,x = step_response(H_highpass)	#Get stepresponse for high pass filter and plot
g1 = General_Plotter("time$\longrightarrow$","$V_o(t)\longrightarrow$","Step response for high pass filter",[],4)
g1.plot_points(t,x)
g1.show()
#********************************************



#Two Frequency Signal
#********************************************
#Plotting input signal
t = np.linspace(0,1e-3,1000000)
g1 = General_Plotter("time$\longrightarrow$","$V_i(t)\longrightarrow$","Two frequency signal",[],5)
g1.plot_points(t,sum_of_sinuosids(t))
g1.show()

#Low pass filter on TFS
t1,y1 = general_input(H_lowpass,sum_of_sinuosids)	#Solve for low pass filter
g1 = General_Plotter("time$\longrightarrow$","$V_o(t)\longrightarrow$","Low pass filter on Two frequency signal",[],6)
g1.plot_points(t1,y1)
g1.show()

#High pass filter on TFS
t1,y1 = general_input(H_highpass,sum_of_sinuosids,1e-5)	#Solve for high pass filter
g1 = General_Plotter("time$\longrightarrow$","$V_o(t)\longrightarrow$","High pass filter on Two frequency signal",[],7)
g1.plot_points(t1,y1)
g1.show()
#********************************************


#Low Frequency damped
#********************************************
#Low pass filter on low freq damped
t = np.linspace(0,0.5,1000000)
t3,y3 = general_input(H_lowpass,damped2,0.5)	#Solve for low pass filter
g1 = General_Plotter("time$\longrightarrow$","Voltage$\longrightarrow$","Low pass filter on low freq damped",[],8)
g1.plot_points(t,damped2(t),"Input signal")	#Plotting input signal
g1.plot_points(t3,y3,"Output signal")		#Plotting output signal
g1.show(1)

#High pass filter on low freq damped
t3,y3 = general_input(H_highpass,damped2,0.5)	#Solve for high pass filter
g1 = General_Plotter("time$\longrightarrow$","Voltage$\longrightarrow$","High pass filter on low freq damped",[],9)
g1.plot_points(t,damped2(t),"Input Signal")	#Plotting input signal
g1.plot_points(t3,y3,"Output Signal")		#Plotting output signal
g1.show(1)
#********************************************


#High frequency damped
#********************************************
#Low pass filter on high freq damped
t = np.linspace(0,1e-3,1000000)
t2,y2 = general_input(H_lowpass,damped1)	#Solve for low pass filter
g1 = General_Plotter("time$\longrightarrow$","damped$\longrightarrow$","Low pass filter on high freq damped",[],10)
g1.plot_points(t,damped1(t),"Input Signal")	#Plotting input signal
g1.plot_points(t2,y2,"Output Signal")		#Plotting output signal
g1.show(1)

#High pass filter on high freq damped
t2,y2 = general_input(H_highpass,damped1)	#Solve for high pass filter
g1 = General_Plotter("time$\longrightarrow$","damped$\longrightarrow$","High pass filter on high freq damped",[],11)
g1.plot_points(t,damped1(t),"Input Signal")	#Plotting input signal
g1.plot_points(t2,y2,"Output Signal")		#Plotting output signal
g1.show(1)
#********************************************

