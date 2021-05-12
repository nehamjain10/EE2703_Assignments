##******************************************************************************************************************************************
#                                                     EE2703 Applied Programming Lab2021
#                                                                Assignment 1

# Purpose  : Reads the netlist file and traverses the circuit definition from last element to first and prints out words in reverse order
# Author   : Neham Jain (EE19B084)
# Input    : Program takes the filename of netlist path as command line input
# Output   : traverses the circuit definition from last element to first and prints out words in reverse order

##******************************************************************************************************************************************

from argparse import ArgumentParser
from sys import exit


#Global Variables declared
circuit_start = ".circuit"
circuit_end = ".end"


def main(params):
    start_index = -1	#Index of line for the .circuit directive
    end_index = -2		#Index of line for the .end directive      
    try:
        with open(params.file_path) as net_list_file:
            lines=net_list_file.readlines()    
            #Finding the values of start_index and end_index, we also check the definition in lower case as someone have might put input in upper case
            for line in lines:
                if str(line[:len(circuit_start)]).lower() == circuit_start:
                    start_index = lines.index(line)
                elif str(line[:len(circuit_end)]).lower() == circuit_end:
                    end_index= lines.index(line)

    except IOError:
        print('No file was found. Please try again')
        exit()
    #checking if file is invalid or not
    if (start_index > end_index or end_index < 0 or start_index < 0):
	    print("Circuit is invalid. Please check the netlist file again.")
	    exit()
    print("Circuit components are: ")
    #Analysing the tokens
    for line in lines[start_index+1:end_index]:	    
        words = line.split('#')[0].split()
        #Determining what type of component is present and printing information
        if words[0][0].lower()=='r':
            print("There is a resistor between nodes {} and {} of value {}".format(words[1],words[2],words[3]))
        if words[0][0].lower()=='l':
            print("There is an resistor between nodes {} and {} of value {}".format(words[1],words[2],words[3]))
        if words[0][0].lower()=='c':
            print("There is a capacitor between nodes {} and {} of value {}".format(words[1],words[2],words[3]))
        if words[0][0].lower()=='v':
            print("There is a voltage source between nodes {} and {} of value {}".format(words[1],words[2],words[3]))
        if words[0][0].lower()=='i':
            print("There is a current source between nodes {} and {} of value {}".format(words[1],words[2],words[3]))


    print("\nPrinting lines of netlist file in reverse order: ")
    #Print the lines in reverse order
    for line in reversed(lines[start_index+1:end_index]):	#Looping the lines in reverse order
	    words = [' '.join(reversed(line.split('#')[0].split()))]	#Removing the comment and then splitting the words in reverse order
	    print(*words)

    
#file is run directly
if __name__=="__main__":
    #Creating an Argument Parser object for parsing command line arguments and providing appropriate help messages
    parser = ArgumentParser(description="Print the reversed lines of the netlist file")
    parser.add_argument("file_path",type=str,help="Provide path for netlist file")
    #Converting parsed arguments into a namespace object and passsing it to the main function
    params=parser.parse_args()
    main(params)
