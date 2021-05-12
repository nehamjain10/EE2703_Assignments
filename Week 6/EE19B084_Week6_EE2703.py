##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 6
#
# Purpose  : Simulating a 1 Dimensional Model of a Tubelight
# Author   : Neham Jain (EE19B084)
# Input    : Values of the grid size, number of electrons, turns, threshold velocity and probability of ionization is given using command line
# Output   : Simulates the electrons in the tubelight, stores the graphs at plots/, saves the tabulated data 

# NOTE: While coding, I have tried to make the code as modular as possible, so that new changes can be introduced without much hassle
##**********************************************************************************************************************************************************************

#Libraries imported in our program
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import pandas as pd
from argparse import ArgumentParser
import os
#Initizialising a random seed, such that every time the user runs the code, we get the same results
np.random.seed(2021)


class General_Plotter():
    ''' Class used for plotting different plots. Shortens the code by quite a bit'''
    
    def __init__(self,xlabel,ylabel,title,fig_num=0):
        ''' xlabel,ylabel,title are used in every graph''' 

        self.xlabel=xlabel
        self.ylabel=ylabel
        self.title=title
        self.fig=plt.figure(fig_num,figsize=(12,10),dpi=100)
        # Increase font size
        matplotlib.rcParams['font.size'] = 16

    
    def general_funcs(self,ax):
        ''' General functions for every graph'''
        
        ax.set_ylabel(self.ylabel)
        ax.set_xlabel(self.xlabel)
        ax.set_title(self.title)
        #self.fig.show()
        self.fig.savefig("plots/"+self.title+".png")

    def population_plot(self,values):
        ''' Plotting histograms'''

        axes=self.fig.add_subplot(111)
        axes.hist(values,bins=100)
        self.general_funcs(axes)

    def scatter_plot(self,X,Y):
        ''' Plotting the phase plots'''

        axes=self.fig.add_subplot(111)
        axes.set_xlim((min(X)-np.std(X),max(X)+np.std(X)))
        axes.set_ylim((min(Y)-np.std(Y),max(Y)+np.std(Y)))
        axes.scatter(X,Y,c='r',marker='X')
        self.general_funcs(axes)



class Simulation():
    ''' Class used for the main simulation . Shortens the code by quite a bit'''
    def __init__(self,params):
        ''' initialising the simulation instance'''

        self.size=params.n
        self.M=params.M
        self.Msig=params.M_sig
        self.timesteps=params.nk
        self.threshold_velocity=params.uo
        self.probability=params.p
        self.dims = self.size*self.M
        self.xx=np.zeros(self.dims)
        self.velocity=np.zeros(self.dims)
        self.displacement=np.zeros(self.dims)
        self.I=[]
        self.X=[]
        self.V=[]

    def calculate_displacement(self,ii):
        ''' Calculating the displacement for each timestep'''

        self.displacement[ii]=self.velocity[ii]+0.5

    def update_attributes(self,ii):
        ''' updating attributes after each timestep'''

        self.calculate_displacement(ii)
        self.xx[ii]=self.xx[ii]+self.displacement[ii]
        self.velocity[ii]=self.velocity[ii]+1
    
    def end_conditions(self,end_particles):
        ''' conditions for when the electron has reached the end of the tubelight'''

        self.xx[end_particles]=0
        self.displacement[end_particles]=0
        self.velocity[end_particles]=0

    def update_displacement_for_collision(self,kl):
        '''Revised model for calculating displacement of collided particles'''  

        time=np.random.rand(len(kl))
        prev_pos=self.xx[kl]-self.displacement[kl]
        displacement=self.velocity[kl]*time+0.5*time*time
        displacement=displacement+0.5*(1-time)**2
        self.xx[kl]=prev_pos+displacement
        self.velocity[kl]=1-time


    def injection(self):
        ''' Injecting electrons into the tubelight at the start of every timestep '''
        m=round(np.random.randn()*self.Msig+self.M)
        free_index=np.where(self.xx==0)[0]
        m=min(m,len(free_index))
        self.xx[free_index[:m]]=1
        return free_index[:m]   

    def update_timestep(self):
        ''' Things to do at every timestep'''

        #Calculate positions where particle is present
        ii=np.where(self.xx>0)[0]
        self.update_attributes(ii)
        
        #Checking particles which have reached the end of tubelight and applying the end conditions
        end_particles=np.where(self.xx>self.size)[0]
        self.end_conditions(end_particles)
        
        #Checking for particles which have crossed threshold velocity
        kk=np.where(self.velocity>=self.threshold_velocity)[0]
        
        #Checking particles that have ionized
        ll=np.where(np.random.rand(len(kk))<=self.probability)[0]
        kl=kk[ll]
        
        #Applying conditions for electrons which have ionized 
        self.update_displacement_for_collision(kl)

        #Storing locations where photon will be emitted
        self.I.extend(self.xx[kl].tolist())

        #Injecting electrons in the tubelight
        free_index=self.injection()
        new_ii=np.concatenate((ii,free_index))

        #Storing positions of the electrons
        self.X.extend(self.xx[new_ii])
        self.V.extend(self.velocity[new_ii])

    def run_simulation(self):
        '''Entry point for running the Simulation '''
        for _ in range(self.timesteps):
            self.update_timestep()
        
        return self.X,self.V,self.I

    def tabulate(self,filename="tabulated_data.txt"):
        '''Tabulating data for intensity vs position'''
        bins = plt.hist(self.I,bins=np.arange(1,self.size,self.size/100))[1]    # Bin positions are obtained
        count = plt.hist(self.I,bins=np.arange(1,self.size,self.size/100))[0]   # Population counts obtained
        
        xpos = 0.5*(bins[0:-1] + bins[1:])     # As no. of end-points of bins would be 1 more than actual no. of bins, the mean of bin end-points are used to get population of count a particular bin
        df = pd.DataFrame()   # A pandas dataframe is initialized to do the tabular plotting of values.
        df['Xpos'] = xpos
        df['count'] = count

        base_filename = filename
        with open(base_filename,'w') as outfile:
            df.to_string(outfile) #Writing the Pandas Dataframe to txt file


def main(params):
    #Creating an instance of the Simulation Class and running the Simulation
    s1=Simulation(params)
    print("Running Simulation...")     
    X,V,I=s1.run_simulation()

    #Plotting the various graphs related to the Simulation
    print("\nPlotting Graphs...")     
    p1=General_Plotter(r"$x$","Number of Electrons","Electron Density",0)
    p1.population_plot(X) 
    p2=General_Plotter(r"$x$","Intensity of Light","Emission Intensity",1)
    p2.population_plot(I)
    p3=General_Plotter(r"$x$","Velocity","Electron Phase Space",2)
    p3.scatter_plot(X,V)
    print("Graphs stored at plots/")     

    #Tabulating Results
    print("\nTabulating Results...")
    s1.tabulate()
    print("Tabulated Files stored at tabulated_data.txt")

#if file is run directly
if __name__=="__main__":
    
    #Creating an Argument Parser object for parsing command line arguments and providing appropriate help messages
    parser = ArgumentParser(description="Specify values of grid size, electrons injected,timesteps, threshold velocity and ionization probability")
    parser.add_argument("-n",type=int,default=100,help="Spatial Grid Size")
    parser.add_argument("-M",type=int,default=5,help="Mean number of electrons injected per turn")
    parser.add_argument("-M_sig",type=int,default=1,help="Std. Deviation of electrons injected per turn")
    parser.add_argument("-nk",type=float,default=500,help="Number of turns to simulate (timesteps)")
    parser.add_argument("-uo",type=float,default=2,help="Threshold Velocity")
    parser.add_argument("-p",type=float,default=0.25,help="Probability that ionization will occur")

    #Converting parsed arguments into a namespace object and passsing it to the main function
    params=parser.parse_args()
    #Creating directory for storing plots
    os.makedirs("plots",exist_ok=True)
    #Running the main function
    main(params)