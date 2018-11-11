#!/usr/bin/python3
import sys

def spin_selection_string(Jmin, Jmax, parities, Nlevels):
	outstring = ""
	for iJ in range(0,int(Jmax-Jmin)+1):
		J = Jmin + iJ
		for pi in parities:
			if J % 1 == 0:
				outstring += "{:.0f}{:s}{:d},".format(J, pi, Nlevels)
			else:
				outstring += "{:.1f}{:s}{:d},".format(J, pi, Nlevels)
	outstring = outstring[0:-1]
	return outstring

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("Usage: ./spin_selection Jmin Jmax parities Nlevels, \ne.g. ./spin_selection 0.5 3.5 +- 10")
		sys.exit(1)
	Jmin = float(sys.argv[1])
	Jmax = float(sys.argv[2])
	if not abs(Jmin % 1 - Jmax % 1) < 0.001:
		raise Exception("Invalid spin range: Jmin =", Jmin, "Jmax =", Jmax)
	parities = str(sys.argv[3])
	if parities == "+":
		parities = ["+"]
	elif parities == "-":
		parities = ["-"]
	elif parities == "+-" or "-+":
		parities = ["+", "-"]
	else:
		raise Exeption("Bad string for parity selection:", parities)
	Nlevels = int(sys.argv[4])

	print(spin_selection_string(Jmin, Jmax, parities, Nlevels))