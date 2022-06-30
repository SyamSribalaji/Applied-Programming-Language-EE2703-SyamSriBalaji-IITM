"""
        EE2703 Applied Programming Lab - 2021
        A1: Spice - Part 1 
        Syam SriBalaji T
        EE20B136
        26/01/22
"""

from sys import argv, exit

"""
To identify beginning and ending of Netlist file by identifying ".circuit" & ".end"
"""
CIRCUIT_START = '.circuit'
CIRCUIT_END = '.end'

"""
To check whether the Python file's name and the Netlist file's name is inputted
"""
if len(argv) != 2:
    print("%s INPUT THE NETLIST FILE" % argv[0])
    exit()

"""
    1. To check whether the given file name is available. If not, it leads to IOError.
    2. If yes, Extracting the circuit definition
    3. Validating the circuit block
    4. Printing the final output

"""
try:
    with open(argv[1]) as NL:
        lines = NL.readlines()
        start = 2; end = 1
        for line in lines:
            if CIRCUIT_START == line[:len(CIRCUIT_START)]:
                start = lines.index(line)
            elif CIRCUIT_END == line[:len(CIRCUIT_END)]:
                end = lines.index(line)
                break
        if start >= end:
            print("Invalid Netlist file")
            exit(0)

        for line in reversed([' '.join(reversed(line.split('#')[0].split())) for line in lines[start+1:end]]):
            print(line)

except IOError:
    print("Invalid Netlist file")
    exit()


"""
Summary of the work done by this code-
    1. Accepting Input of Netlist file
    2. Verifying whether the file is available
    3. Checking whether file has correct format
    4. Removing the comments 
    5. Print the reversed form of the text in Netlist file as Output
"""