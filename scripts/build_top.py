#!/usr/bin/python3
import os
import argparse

parser = argparse.ArgumentParser(description = "Take a .gro file, reorder the file such that all molecules are together and in order, then cut out any extra unwanted molecules specified with -l.\
For example, <DOPC:5> will remove 5 entire DOPC molecules.",
 epilog = "Using the optional arguments -p and -po you may rewrite a topology file automatically. ")

parser.add_argument("-c", "--input", help = "<.gro> input structure file", required = True)
parser.add_argument("-p", "--topology", help = "<.top> input topology file")
parser.add_argument("-po", "--top_out", help = "<.top> output topology file", default = "new_topol.top")

args = parser.parse_args()

# get the names of all molecules in the gro file
mol_list = []
with open(args.input, 'r') as f:
    for line in f.readlines()[2:-1]:
        if line[5:11] not in mol_list:
            mol_list.append(str(line[5:11]))

print(mol_list)

# first need first index of each molecule:
first_indices = []
for mol in mol_list:
    with open(args.input, "r") as f:
        for line in f.readlines()[2:-1]:
            #capture first instance of numbered mol
            if mol in line:
                first_indices.append(str(line[0:11]))
                break
# number of beads
beads_list = []
for index in first_indices:
    counter = 0
    with open(args.input, "r") as f:
        for line in f.readlines()[2:-1]:
            if index in line:
                counter += 1
    beads_list.append(counter)

beads_dict = {}
for i in range(len(beads_list)):
    beads_dict[mol_list[i].rstrip()] = beads_list[i]

if args.topology is not None:
    writer = open(args.top_out, "w")
    with open(args.topology, "r") as f:
        for line in f.readlines():
            if "molecules" not in line:
                writer.write(line)
            else:
                break
    writer.write("[ molecules ]\n; name\tnumber\n")

for mol in mol_list:
    counter = 0
    with open(args.input, "r") as f:
        for line in f.readlines()[2:-1]:
            if mol in line:
                counter += 1
    if mol.rstrip() == "ION":
        # print("NA\t\t" + str(int(counter / beads_dict[mol.rstrip()])) + "\n")
        writer.write("NA\t\t" + str(int(counter / beads_dict[mol.rstrip()])) + "\n")
    else:
        print(mol.rstrip() + "\t\t" + str(int(counter / beads_dict[mol.rstrip()])) + "\n")
        writer.write(mol.rstrip() + "\t\t" + str(int(counter / beads_dict[mol.rstrip()])) + "\n")