##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 5
#
# Purpose  : Solving, Laplace's Equation for values of the potential and current
# Author   : Neham Jain (EE19B084)
# Input    : Values of nx,ny,radius of wire and iterations to perform can be specified using the command line
# Output   : Plots different graphs as specified in the report
#
# NOTE: While plotting, we have tried to reduce the code redundancy as much as possible to minimize the size of the code
##**********************************************************************************************************************************************************************

#Libraries imported in our program
import numpy as np 
import matplotlib.pyplot as plt 
from argparse import ArgumentParser
import mpl_toolkits.mplot3d.axes3d as p3
import os

class General_Plotter():
    ''' Class used for plotting different plots. Shortens the code by quite a bit'''
    
    def __init__(self,xlabel,ylabel,title,zlabel=None):
        ''' xlabel,ylabel,title are used in every graph''' 

        self.xlabel=xlabel
        self.ylabel=ylabel
        self.zlabel=zlabel
        self.title=title
        self.fig=plt.figure()

    def general_funcs(self,ax):
        ''' General functions for every graph'''
        
        ax.set_ylabel(self.ylabel)
        ax.set_xlabel(self.xlabel)
        ax.set_title(self.title)
        self.fig.savefig("plots/"+self.title+".png")

    '''Specialized functions for plotting different types of graphs'''

    def contour_plot(self,x,y,phi,wire_loc,cmap=plt.get_cmap('hot')):
        axes=self.fig.add_subplot(111)
        t=axes.contourf(y,x[::-1],phi,cmap=cmap)
        axes.plot(x[:,0][wire_loc[0]],y[0,:][wire_loc[1]],'ro')
        axes.grid()
        self.fig.colorbar(t)
        self.general_funcs(axes)

    def general_plot(self,x1,y1):
        axes=self.fig.add_subplot(111)
        axes.plot(x1,y1,'-r',markersize=3,label='original')
        self.general_funcs(axes)

    def semilogy_plot(self,x1,y1,x2=None,y2=None,x3=None,y3=None):
        axes=self.fig.add_subplot(111)
        axes.semilogy(x1,y1,'ro',label='Original')
        if x2 is not None:
            axes.semilogy(x2,y2,'g-',label='Fit1')
            axes.semilogy(x3,y3,'b-',label='Fit2 (>500 iterations)')
            axes.legend()
        self.general_funcs(axes)
    
    def loglog(self,x1,y1,x2=None,y2=None,x3=None,y3=None):
        axes=self.fig.add_subplot(111)
        axes.loglog(x1,y1,'ro',markersize=3,label='Original')
        if x2 is not None:
            axes.loglog(x2,y2,'g-',markersize=3,label='Fit1')
            axes.loglog(x3,y3,'b-',markersize=3,label='Fit2 (>500 iterations)')
        axes.legend()
        self.general_funcs(axes)
    
    def threeD_plot(self,X,Y,phi,cmap=plt.get_cmap('jet')):
        axes=p3.Axes3D(self.fig) 
        surf = axes.plot_surface(Y, X, phi.T, rstride=1, cstride=1,cmap=cmap)
        self.fig.colorbar(surf,shrink=0.5,pad=0.08)
        axes.set_zlabel(self.zlabel,fontsize=15) 
        self.general_funcs(axes)

    def quiver_plots(self,X,Y,Jx,Jy,wire_loc):
        axes=self.fig.add_subplot(111)
        axes.quiver(Y,X[::-1],Jx,Jy, scale = 5)
        axes.plot(X[:,0][wire_loc[0]],Y[0,:][wire_loc[1]],'ro')
        self.general_funcs(axes)
    

def create_grid(x,y,radius):    
    ''' Creating grid of coordinates of the plate and initialising phi'''
    phi_initial=np.zeros((y,x))
    xx = np.linspace(-0.5,0.5,y)
    yy = np.linspace(-0.5,0.5,x)
    Y_coords,X_coords = np.meshgrid(yy,xx)
    wire_loc = np.where((X_coords**2 + Y_coords**2) < (radius/x)**2)
    phi_initial[wire_loc]=1
    return X_coords,Y_coords,phi_initial,wire_loc

def Potential_Solver(phi,epochs,wire_loc):
    '''Function to solve potentials at various x,y locations in the plates'''
    
    errors = np.zeros(epochs)
    def boundary_conditions(phi,wire_loc):
        ''' Function to enforce appropriate boundary conditions''' 
    
        phi[1:-1,0] = phi[1:-1,1]   #Left side boundary condition
        phi[0,1:-1] = phi[1,1:-1]   #Top side
        phi[1:-1,-1] = phi[1:-1,-2] #Right side boundary condition
        phi[-1,1:-1] = 0            #Bottom part of plate is kept grounded
        phi[wire_loc] = 1.0         #Wire locations are at 1V
        return phi

    for i in range(epochs):
        oldphi=phi.copy()   #Need a proper copy not a view
        phi[1:-1,1:-1]=0.25*(phi[1:-1,0:-2]+ phi[1:-1,2:]+phi[0:-2,1:-1]+ phi[2:,1:-1])    #Laplace Equation
        phi=boundary_conditions(phi,wire_loc)
        errors[i]=(abs(phi-oldphi)).max()   #Calculate errors.
    
    return phi,errors


def solve_currents(phi,Nx,Ny):
    '''
    Function computes the current vectors ie: Jx , Jy
    '''
    Jx = np.zeros((Ny,Nx))
    Jy = np.zeros((Ny,Nx))
    Jx[:,1:-1] = 0.5*(phi[:,0:-2]-phi[:,2:])
    Jy[1:-1,:] = 0.5*(phi[2:, :]-phi[0:-2,:])

    return Jx,Jy

def solve_temp(Ny,Nx,Jx,Jy,iters,wire_loc):
    '''
    Function computes the temperature vectors
    '''
    T=300*np.ones((Ny,Nx))
    for i in range(iters    ):
        T[1:-1,1:-1]=0.25*(T[1:-1,0:-2]+ T[1:-1,2:]+ T[0:-2,1:-1] + T[2:,1:-1]+(Jx[1:-1,1:-1])**2 +(Jy[1:-1,1:-1])**2)
        T[1:-1,0]=T[1:-1,1] # left boundary
        T[1:-1,-1]=T[1:-1,-2] # right boundary
        T[0,1:-1]=T[1,1:-1] # top boundary
        T[-1,1:-1] =300 # ground
        T[wire_loc]=300.0 #wire is at 300K
    return T

def fit_exp(x,A,B):
    ''' Function to obtain an exponential '''
    return A*np.exp(B*x)

def error_fit(x,y):
    ''' Evaluates the parameters of an exponent. x is the data vector, y is the vector of true values ''' 
    logy=np.log(y)
    xvec=np.zeros((len(x),2))
    xvec[:,0]=x
    xvec[:,1]=1
    B,logA=np.linalg.lstsq(xvec, np.transpose(logy),rcond=None)[0]
    return (np.exp(logA),B)

def max_error(A,B,N):
    ''' Finds an upper bound for the error ''' 
    return -A*(np.exp(B*(N+0.5)))/B


def main(params):
    
    #Parsing arguments
    x_length=params.x
    y_length=params.y
    radius=params.r
    iters=params.iters
    
    #Creating the mesh grid of the plate
    X_coords,Y_coords,phi_initial,wire_loc=create_grid(x_length,y_length,radius)
    
    #Plotting contour plot of the initial potential assumption of the plate
    p1=General_Plotter(xlabel="X-axis",ylabel="Y-axis",title="Contour Plot of Initial Potential Assumption")
    p1.contour_plot(X_coords,Y_coords,phi_initial,wire_loc)
    
    #Solving the potential
    print("Solving Potential.... ")
    phi,error=Potential_Solver(phi_initial,iters,wire_loc)
    
    #Fitting an exponential to the error data
    print("Plotting Error Values.... ")
    A,B = error_fit(range(iters),error) # fit1
    A_500,B_500 = error_fit(range(iters)[500:],error[500:]) # fit2
    
    #Plotting the evolution of the error function 
    p2=General_Plotter(xlabel=r'iters$\rightarrow$',ylabel=r'Error$\rightarrow$',title='Plot of Error vs number of iterations')
    p2.general_plot(range(iters),error)
    
    #Error, fit1 and fit2 in a semilog plot 
    p3=General_Plotter(xlabel=r'iters$\rightarrow$',ylabel=r'Error$\rightarrow$',title='Semilog plot of Error vs number of iterations')
    p3.semilogy_plot(range(iters)[::50],error[::50],range(iters)[::50],fit_exp(range(iters)[::50],A,B),range(iters)[::50],fit_exp(range(iters)[::50],A_500,B_500))
    
    #Error, fit1 and fit2 in a loglog plot
    p4=General_Plotter(xlabel=r'Niters$\rightarrow$',ylabel=r'Error$\rightarrow$',title='Loglog plot of Error vs number of iterations')
    p4.loglog(range(iters)[::50],error[::50],range(iters)[::50],fit_exp(range(iters)[::50],A,B),range(iters)[::50],fit_exp(range(iters)[::50],A_500,B_500))

    #Plotting max error in semilog plot  
    p5=General_Plotter(xlabel=r'Niters$\rightarrow$',ylabel=r'Error$\rightarrow$',title='Semilog plot of Cumulative Error vs number of iterations')
    p5.semilogy_plot(range(iters)[::50],max_error(A,B,np.arange(0,iters,50)))
 
    #3-D plot of potential
    p6=General_Plotter(xlabel=r'x$\rightarrow$',ylabel=r'y$\rightarrow$',title='The 3-D surface plot of the potential',zlabel=r'$\phi\rightarrow$')
    p6.threeD_plot(X_coords,Y_coords,phi)

    #Contour plot of the final potential
    p7=General_Plotter(xlabel="X-axis",ylabel="Y-axis",title="Contour Plot of Actual Potential")
    p7.contour_plot(X_coords,Y_coords,phi,wire_loc)

    #Solving to obtain the currents Jx and Jy
    print("Solving Current.... ")
    Jx,Jy=solve_currents(phi,x_length,y_length)
    
    #Quiver plot of the current
    p8=General_Plotter(xlabel=r"x$\rightarrow$",ylabel=r'y$\rightarrow$',title="Current Density")
    p8.quiver_plots(X_coords,Y_coords,Jx,Jy,wire_loc)

    #Solving for temperature in the plates
    print("Solving Temperature.... ")
    T=solve_temp(y_length,x_length,Jx,Jy,iters,wire_loc)

    #Plotting the 3D Plot for temperature
    p9=General_Plotter(xlabel=r'x$\rightarrow$',ylabel=r'y$\rightarrow$',title='The 3-D surface plot of the Temperature',zlabel=r'$Temperature \rightarrow$')
    p9.threeD_plot(X_coords,Y_coords, T)
    print("Plots saved at plots/")    

#if file is run directly
if __name__=="__main__":
    
    #Creating an Argument Parser object for parsing command line arguments and providing appropriate help messages
    parser = ArgumentParser(description="Specify values of length along x, length along y, radius of central lead, number of iters. to perform")
    parser.add_argument("-x",type=int,default=25,help="Length of plate along x")
    parser.add_argument("-y",type=int,default=25,help="Length of plate along y")
    parser.add_argument("-r",type=float,default=8,help="Radius of central lead")
    parser.add_argument("-iters",type=int,default=1500,help="Number of iterations to perform")

    #Converting parsed arguments into a namespace object and passsing it to the main function
    params=parser.parse_args()
    #Creating directory for storing plots
    os.makedirs("plots",exist_ok=True)
    #Running the main function
    main(params)
