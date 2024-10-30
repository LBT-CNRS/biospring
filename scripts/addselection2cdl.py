#!/usr/bin/python

import os
import sys

if len(sys.argv) < 4:
    print "Usage: addselection2cdl.py cdlinput cdloutput sel1name [sel2name sel3name ...]"

#Read selection names and ask for atom ids
selections = {}
for n in sys.argv[3:]:
    name = n + "_selection"
    ids = raw_input("Atom ids for selection %s? (separate by space)\n"% name)
    selections[name] = [int(x) for x in ids.strip().split()]

#Copy cdlin file and adds the selection fields    
cdlin = open(sys.argv[1],"r")
cdlout = open(sys.argv[2],"w")

#dimensions field
line = cdlin.readline()
while line.find("variables") == -1:
    cdlout.write(line)
    line = cdlin.readline()
    
for s in selections.keys():
    cdlout.write("\t%s = %i;\n" % (s, len(selections[s])))
cdlout.write("\n")

#variables field
while line.find("data:") == -1:
    cdlout.write(line)
    line = cdlin.readline()

for s in selections.keys():
    cdlout.write("\tint	%s(%s);\n" %  (s,s))
    cdlout.write("\t\t%s:long_name = \"Selection particle ids list\";\n\n" % s)
cdlout.write("\n")

#data field
while line.strip() != "}":
    cdlout.write(line)
    line = cdlin.readline()

for s in selections.keys():
    cdlout.write("\t%s = \n" % s)
    for i in selections[s][:-1]:
        cdlout.write("%i,\n" % i)
    cdlout.write("%i;\n" % selections[s][-1])
    
while line:
    cdlout.write(line)
    line = cdlin.readline()
    
cdlout.close()
cdlin.close()

for s in selections.keys():
    print "added selection %s %s" % (s, selections[s])


