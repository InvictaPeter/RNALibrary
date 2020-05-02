import time

import sys, os, argparse, subprocess, re

def stdout_from_command(command):
    p = subprocess.Popen(command, stdout = subprocess.PIPE, shell = True)
    return iter(p.stdout.readline, b'')

f = open("FinalBatch.txt", "r")

def file_len(fname):
    with open(fname) as g:
        for i, l in enumerate(g):
            pass
    return i + 1
z=[]
for x in range(64224*3+2):
	temp=f.readline()
	if(temp[0]=="."):
		f.readline()
		temp2=f.readline().split()
		z.append(temp2[0])

def get_rna_structure_energy(sequence): #From Original Program
    if len(sequence) < 1:
        return (0, 0)

    lines = stdout_from_command("echo %s | RNAfold --noPS --MEA -p" % sequence)
    
    #first line is just the seq - skip
    lines.next()

    #second line has the MFE
    MFE = float(re.sub('[()]','', lines.next().split()[-1]))

    #third line has the ensemble energy - discarding for now
    lines.next()

    #fourth line has the centroid energy
    nextline = lines.next().split()
    nrg = re.sub('[{}]','', nextline[-2])
    strx = nextline[0]

    #fifth line has the MEA
    MEA = float(re.sub('[{}]','', lines.next().split()[-2]))

    return (nrg)
def get_moststable(sequence): #RNALfold Stable Structure Finder
    a=[]
    if len(sequence) < 1:
        return (0, 0)
    lines = stdout_from_command("echo %s | RNALfold --span=50" % sequence)
    for i in range(200): #Get all the outputs. Probably not over 200. If it is, make this number bigger rinse repeat
        try:
            p=lines.next()
            a.append(p)
        except:
            pass
    stableseq=(a[-2])[0:-1]
    stableenergy=(float((a[-1])[2:-2]))
    return stableseq, stableenergy
c=[]
for b in z:
    if (len(b)>150):
        c.append(b)
# length (which should be constant)
# folding free energy for the entire 5' UTR
# folding free energy for the most stable 50 nucleotide stretch (which can be computed using another ViennaRNA program called RNALfold)
# number of upstream NUG codons
# number of upstream stop codons
# GC content (percent of the UTR that is G/C nucleotides)
def countNUG(inseq):
	return inseq.count('AUG')+inseq.count('UUG')+inseq.count('CUG')+inseq.count('GUG')
def findStop(inseq):
	return inseq.count("UAA")+inseq.count("UGA")+inseq.count("UAG")
def GCcontent(inseq):
	return float((inseq.count("G")+inseq.count("C")))/len(inseq)*100
def GenMetrics(inseq): 
	print("5' UTR length: "+ str(len(inseq)))
	print("5' UTR folding free energy: "+str(get_rna_structure_energy(inseq)))
	print("Lowest MFE of sequence of length 50: " + str(get_moststable(inseq)[1]))
	print("Amount of upstream NUG occurances: "+str(countNUG(inseq)))
	print("Amount of upstream stop codon occurances: "+str(findStop(inseq)))
	print("Percent GC content: "+str(GCcontent(inseq)))
GenMetrics(c[23]) 


