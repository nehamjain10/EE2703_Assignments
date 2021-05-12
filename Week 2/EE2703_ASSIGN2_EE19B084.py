##**********************************************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab 2021
#                                                                 Assignment 2

# Purpose  : Reads the netlist file, traverses the circuit definition and solves for the Nodal Voltages
# Author   : Neham Jain (EE19B084)
# Input    : Program takes the filename of netlist path as command line input
# Output   : Prints all the unknowns (the nodal voltages and the current through the voltage sources)

##**********************************************************************************************************************************************************************

'''
------BASIC SYNTAX OF THE NETLIST FILE------
.circuit	#This is a comment
V1 node1 node2 dc value 			#DC VOLTAGE SOURCE -- (node1 == +ve terminal) (node2 == -ve terminal)
V2 node1 node2 ac Vp-p phase 		#AC VOLTAGE SOURCE -- Vp-p is peak to peak voltage, phase in degrees
R1 node1 node2 value 				#RESISTANCE
L1 node1 node2 value 				#INDUCTOR
C1 node1 node2 value 				#CAPACITOR
I1 node1 node2 value 				#CURRENT SOURCE -- arrow from node1 to node2
.end
.ac V__ 						#frequency of AC voltage source -- Currently only one AC voltage source is supported


Also the GND node is compulsory. It has to be in CAPITAL LETTERS.
While printing the output, it is rounded to 5 decimal places
'''

#Importing necessary libraries to be used in our program
from argparse import ArgumentParser
from sys import exit
import numpy as np
import math
import cmath
import warnings
warnings.filterwarnings('ignore')

#Global Variables to be used in the program
CIRCUIT_START = ".circuit"
CIRCUIT_END = ".end"
AC=".ac "
ac_flag=False      

#List of nodes,passive elements,independent sources and dependent courses
node_list=[]
passive_element_list=[]
independent_src_list=[]
matrix_size=0

#Class for all passive elements and storing their relevant parameters
class Passive_Element():
    global ac_flag
    def __init__(self,name,type,value,node1,node2,w):
        #Identifier of that particular element
        self.name=name
        #Storing impedance of that linear elements such as resistor, capacitor, inductor
        if type=="r":
            self.impedance=value
        if type=='c':
            if ac_flag==True:
                self.impedance=1/(1j*w*value)
            else:
                self.impedance=1e15
        if type=='l':
            if ac_flag==True:
                self.impedance=1j*w*value
            else:
                self.impedance=1e-15

        #Storing the nodes that the particular element is connnected to
        self.node1=node1
        self.node2=node2

#Class for all storing all nodes and connected elements to that particular node
class IndependentSrc():
    def __init__(self,name,type,node1,node2,ac_or_dc,value,phase=0):
        global matrix_size
        self.name=name
        self.type=type
        #Current through that particular voltage is an unknown so we increase the size of matrix
        if type=='v':
            matrix_size+=1
        self.node1=node1
        self.node2=node2
        self.ac_or_dc=ac_or_dc
        if ac_flag==True:
            self.value=value*cmath.exp(1j*phase)/2
        else:
            self.value=value

#Class for all storing all nodes and connected elements to that particular node
class Node():
    def __init__(self,name):
        global matrix_size
        self.passive_connections={nodes:[] for nodes in node_list}
        self.sources={nodes:[] for nodes in node_list}
        self.name=name
        self.nodal_voltage=0
        matrix_size+=1
        for i in passive_element_list:
            if i.node1==name:
                self.passive_connections[i.node2].append(i.impedance)
            elif i.node2==name:
                self.passive_connections[i.node1].append(i.impedance)
        
        for i in independent_src_list:
            if i.node1==name:
                self.sources[i.node2].append((i.name,i.type,i.value,'p'))
            elif i.node2==name:
                self.sources[i.node1].append((i.name,i.type,i.value,'n'))


#Function to find the relevant portion and discards comments and other garbage part of the circuit
def find_circuit_portion(params):
    start_index = -1	#Index of line for the .circuit directive
    end_index = -2		#Index of line for the .end directive
    frequency=0
    global ac_flag
    try:
        with open(params.file_path) as net_list_file:
            lines=net_list_file.readlines()    
            #Finding the values of start_index and end_index, we also check the definition in lower case as someone have might put input in upper case
            for line in lines:
                if str(line[:len(CIRCUIT_START)]).lower() == CIRCUIT_START:
                    start_index = lines.index(line)
                elif str(line[:len(CIRCUIT_END)]).lower() == CIRCUIT_END:
                    end_index= lines.index(line)
                elif str(line[:len(AC)]).lower() == AC:
                    ac_flag = True
                    frequency=int(line.split()[2])
    except IOError:
        print('No file was found. Please try again')
        exit()
    #checking if file is invalid or not
    if (start_index > end_index or end_index < 0 or start_index < 0):
	    print("Circuit is invalid. Please check the netlist file again.")
	    exit()
    return lines,start_index,end_index,frequency


#Parsing the relevant portion in the circuit and storing all the components of the circuit
def parse_file(lines,start_index,end_index,frequency):
    global ac_flag
    if ac_flag==True:
        w=2*math.pi*frequency
    else:
        w=0
    #parsing the lines and creating the relevant objects of the class based on the type of element
    for line in (lines[start_index+1:end_index]):	#Looping the lines in reverse order
        line = line.split('#')[0]	#Removing the comment
        components=line.split()
        type=components[0][0].lower()
        name=components[0]
        #If type of element is Passive
        if type=='r' or type=='l' or type=='c':
            node1=components[1]
            node2=components[2]
            value=float(components[3])
            #storing Passive_Element object in a list containing all the passive objects
            passive_element_list.append(Passive_Element(name,type,value,node1,node2,w))
            if node1 not in node_list:
                node_list.append(node1)
            if node2 not in node_list:
                node_list.append(node2)
        #If the type of element is Independent Source, creating the object of type IndependentSrc
        if type=='v' or type=='i':
            node1=components[1]
            node2=components[2]
            ac_or_dc=components[3]
            value=float(components[4])
            phase=0
            if node1 not in node_list:
                node_list.append(node1)
            if node2 not in node_list:
                node_list.append(node2)
            if ac_or_dc=='ac':
                phase=float(components[5])
            independent_src_list.append(IndependentSrc(name,type,node1,node2,ac_or_dc,value,phase))   

#Function to make equations regarding the modified nodal analysis
def make_equations(nodal_details):
    #creating a matrix containing complex values or float depending on whether circuit is ac or dc
    if ac_flag==True:
        A=np.zeros((matrix_size,matrix_size),dtype=np.complex128)
        b=np.zeros(matrix_size,dtype=np.complex128)
    else:
        A=np.zeros((matrix_size,matrix_size),dtype=np.float64)
        b=np.zeros(matrix_size,dtype=np.float64)
    independent_vol_src=[i for i in independent_src_list if i.type=='v']
    nodes=list(nodal_details.keys())
    unknowns=nodes.copy()
    names_vol_src=[i.name for i in independent_vol_src]
    unknowns.extend(names_vol_src)
    #Creating a list called unknowns which contains values of all the unknown parameters whose value we have to determine
    count=0

    #Parsing through all the nodes of the circuits and performing modified nodal analysis to form our matrix
    for i in nodes:
        coeffs={uk:0 for uk in unknowns}
        passive=nodal_details[i].passive_connections
        if i=="GND":
            coeffs[i]=1
            b[count]=0
            A[count,:]=list(coeffs.values())
            count+=1
            continue
        for vals in passive.items():
            if len(vals[1])==0:
                continue
            temp=0
            for k in vals[1]:
                temp=temp+np.float64(1)/k
            coeffs[vals[0]]=-1*temp
            coeffs[i]=coeffs[i]+temp
        src_list=nodal_details[i].sources
        for src in src_list.items():
            if len(src[1])>0:
                for source in src[1]:
                    if source[1]=='v':
                        coeffs[source[0]]=1 if source[3]=='n' else -1                        
                    if source[1]=='i':
                        multiplier = 1 if source[3] == 'n' else -1
                        b[count]=b[count]+multiplier*source[2]
        A[count,:]=list(coeffs.values())
        count+=1

    #Parsing through all the voltage sources to add equations corresponding to that in our matrix
    for src in independent_vol_src:
        coeffs={uk:0 for uk in unknowns}
        coeffs[src.node1]=1
        coeffs[src.node2]=-1
        A[count,:]=list(coeffs.values())
        b[count]=src.value
        count+=1
    return A,b,names_vol_src

#function to print all the unknowns in the circuit and also writing the output to a file named output.txt
def print_results(names_vol_src,node_list,sols,path):
    file1 = open("output.txt","w")  
    count=0
    print("\n")
    for i in node_list:
        temp=f"The voltage at node {i} is: {round(sols[count],5)} V"
        file1.write(temp+"\n") 
        print(temp)

        count+=1
    for i in names_vol_src:
        temp=f"The current through voltage source {i} is: {round(sols[count],5)} A"
        file1.write(temp+"\n") 
        print(temp)
        count+=1
    file1.close()

#main function to be run
def main(params):
    #Finding the relevant function
    lines,start_index,end_index,frequency=find_circuit_portion(params)
    
    #Parsing the function
    parse_file(lines,start_index,end_index,frequency)

    if "GND" not in node_list:
        print("No ground node specified in netlist file. Please check and try again.")
        exit()
    
    nodal_details={i:None for i in node_list}
    
    for i in node_list:
        nodal_details[i]=Node(i)
    
    #Making equations from the nodal details
    A,b,names_vol_src=make_equations(nodal_details)
    
    #using solve function of linalg library to solve our matrix
    sols=np.linalg.solve(A,b)
    
    if ac_flag==False:
        sols=sols.astype(np.float64)
    print_results(names_vol_src,node_list,sols,params.file_path)


#file is run directly
if __name__=="__main__":
    #Creating an Argument Parser object for parsing command line arguments and providing appropriate help messages
    parser = ArgumentParser(description="Print the reversed lines of the netlist file")
    parser.add_argument("file_path",type=str,help="Provide path for netlist file")
    #Converting parsed arguments into a namespace object and passsing it to the main function
    params=parser.parse_args()
    main(params)