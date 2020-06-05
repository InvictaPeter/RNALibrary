import time

import sys, os, argparse, subprocess, re
z=[]
count=[]
metriclist=[]
def stdout_from_command(command):
    p = subprocess.Popen(command, stdout = subprocess.PIPE, shell = True)
    return iter(p.stdout.readline, b'')

def file_len(fname):
    with open(fname) as g:
        for i, l in enumerate(g):
            pass
    return i + 1


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

def countNUG(inseq):
	return inseq.count('AUG')+inseq.count('UUG')+inseq.count('CUG')+inseq.count('GUG')
def findStop(inseq):
	return inseq.count("UAA")+inseq.count("UGA")+inseq.count("UAG")
def GCcontent(inseq):
	return float((inseq.count("G")+inseq.count("C")))/len(inseq)*100
def verifystructure(sequence,structure): #From Original Program
    if len(sequence) < 1:
        return (0, 0)
    lines = stdout_from_command("echo %s | RNAfold --noPS --MEA -p" % sequence)
   
    #first line is just the seq - skip
    lines.next()
    secstruct=lines.next().split()[0]
    if(secstruct==structure):
        return 'Verified'
    else:
        return 'Discard'
def GenMetrics(inseq,instruct): 
	metrics=[]
	metrics.append(instruct)
	metrics.append(inseq)
	metrics.append("5' UTR length: "+ str(len(inseq)))
	metrics.append("5' UTR folding free energy: "+str(get_rna_structure_energy(inseq)))
	metrics.append("Lowest MFE of sequence of length 50: " + str(get_moststable(inseq)[1]))
	metrics.append("Amount of upstream NUG occurances: "+str(countNUG(inseq)))
	metrics.append("Amount of upstream stop codon occurances: "+str(findStop(inseq)))
	metrics.append("Percent GC content: "+str(GCcontent(inseq)))
	metrics.append("Match: "+ verifystructure(inseq,instruct))
	
	return metrics
def RetrieveFile(halffilelen):
	f = open("testrun.txt", "r")
	
	counter=0
	for x in range(halffilelen*2):
		temp=f.readline()
		if(temp[0]=="."):
			structure=temp.split()
			z.append(structure[0])
			count.append(counter)
			f.readline()
			f.readline()
			targetseq=f.readline().split()#final line of the rnainverse command is usually best sequence
			z.append(targetseq[0])
		counter+=1	
	targetstructs=z[0::2]
	counter=0
	for item in z[1::2]:
		metriclist.append(GenMetrics(item,targetstructs[counter]))
		counter+=1
		
	f.close()
def GenMetricFile(diversity):
	f=open("Metrics.txt", "w")
	iterable=0
	for item in metriclist:
		for item2 in item:
			print(item2)
			f.write(item2)
			f.write('\n')
		f.write("diversity: "+str(diversity[iterable]))
		iterable+=1
		f.write('\n')
		f.write('\n')
		
			



