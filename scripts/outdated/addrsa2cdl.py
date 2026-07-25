#!/usr/bin/python

import sys
import os

cdl = open(sys.argv[1],"r")
rsa = open(sys.argv[2],"r")
out = open(sys.argv[1] + ".out","w")

#Go to the surfaceaccessibility field in the cdl
lineCDL = cdl.readline()
while lineCDL.find("surfaceaccessibility =") == -1 :
	out.write(lineCDL)
	lineCDL = cdl.readline()


#Go to the beginning of the RES SAS table	
line = rsa.readline().split()
while line[0] != "RES":
	line = rsa.readline().split()

#Copy the ABS to the cdl
#First residue	
sasAbs = line[5]
out.write( lineCDL + sasAbs )
line = rsa.readline().split()
while line[0] == "RES":
	sasAbs = line[5]
	out.write(",\n" + sasAbs)
	line = rsa.readline().split()
#Finish with ;
out.write(";\n")

#Go to the next field in the cdl
lineCDL = cdl.readline()
while lineCDL.find(";") == -1 :
	lineCDL = cdl.readline()

while lineCDL :
	lineCDL = cdl.readline()
	out.write(lineCDL)


rsa.close()
cdl.close()
out.close()

os.remove(sys.argv[1])
os.rename(sys.argv[1] + ".out",sys.argv[1])
