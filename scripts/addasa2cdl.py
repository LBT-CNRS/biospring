#!/usr/bin/python

import sys
import os
import argparse


def addASA(args):

	cdl = open(args.cdl,"r")
	asa = open(args.asa,"r")
	out = open(args.cdl + ".out","w")

	sasaCol=-1
	if(args.prog=="freesasa"):
		sasaCol=9
	elif(args.prog=="naccess"):
		sasaCol=8
	else:
		print("error column sasa")
		quit()

	# Go to the surfaceaccessibility line in the cdl
	# while writing to output
	lineCDL = ""
	while lineCDL.find("surfaceaccessibility =") == -1 :
		lineCDL = cdl.readline()
		out.write(lineCDL)

	# Copy the ABS to the cdl
	sasas = [line.strip().split()[sasaCol] for line in asa.readlines() if line.startswith("ATOM")]
	for s in sasas[:-1]:
		out.write(s + ",\n")
	out.write(sasas[-1] + ";\n")

	# Go to the next field in the cdl
	lineCDL = cdl.readline()
	while lineCDL.find(";") == -1 :
		lineCDL = cdl.readline()

	# Write the rest
	while lineCDL :
		lineCDL = cdl.readline()
		out.write(lineCDL)

	asa.close()
	cdl.close()
	out.close()

	os.remove(args.cdl)
	os.rename(args.cdl + ".out", args.cdl)


parser = argparse.ArgumentParser(description="Add ASA values to CDL file.")
parser.add_argument('cdl', type=str, help="CDL input file")
parser.add_argument('asa', type=str, help="ASA input file")
parser.add_argument('-p','--prog', type=str, choices=['naccess', 'freesasa'], help="Program used to generate ASA file.")

if __name__ == '__main__':
	args = parser.parse_args(sys.argv[1:])
	print(args.cdl)
	addASA(args)