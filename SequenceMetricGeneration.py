import time

import sys, os, argparse, subprocess, re, Levenshtein
structure_and_sequence_list=[]
seqnumber=0
metriclist=[]
EditDistanceList=[]
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
	lines.__next__()

	#second line has the MFE
	MFE = float(re.sub('[()]','', lines.__next__().decode("utf-8").split()[-1]))

	#third line has the ensemble energy - discarding for now
	lines.__next__()

	#fourth line has the centroid energy
	nextline = lines.__next__().decode("utf-8").split()
	nrg = re.sub('[{}]','', nextline[-2])
	strx = nextline[0]

	#fifth line has the MEA
	MEA = float(re.sub('[{}]','', lines.__next__().decode("utf-8").split()[-2]))

	return (nrg)
def get_moststable(sequence): #RNALfold Stable Structure Finder
	a=[]
	if len(sequence) < 1:
		return (0, 0)
	lines = stdout_from_command("echo %s | RNALfold --span=50" % sequence)
	for i in range(seqnumber*2): #Get all the outputs. Probably not over 200. If it is, make this number bigger rinse repeat
		#print(a)
		try:
			a.append(lines.__next__())
			#print(a)
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
	lines.__next__()
	secstruct=lines.__next__().decode("utf-8").split()[0]
	if(secstruct==structure):
		return 'Verified'
	else:
		return 'Discard'

def Levenshtein_distance(structures):
	distances=[]
	for i in range(1,len(structures)+1):
		candidate=structures[i-1]
		restofseq=[structures[0:i-1]+structures[i:]]
		diversitysum=0
		for item in restofseq[0]:
			diversitysum+=Levenshtein.distance(candidate,item)
		distances.append(float(diversitysum/float(str(len(structures))+".0")))
	return distances

def GenMetrics(inseq,instruct,counter):
	metrics=[]
	metrics.append(instruct)
	metrics.append(inseq)
	metrics.append("5' UTR length: "+ str(len(inseq)))
	metrics.append("5' UTR folding free energy: "+str(get_rna_structure_energy(inseq)))
	metrics.append("Lowest MFE of sequence: " + str(get_moststable(inseq)[1]))
	metrics.append("Amount of upstream NUG occurances: "+str(countNUG(inseq)))
	metrics.append("Amount of upstream stop codon occurances: "+str(findStop(inseq)))
	metrics.append("Percent GC content: "+str(GCcontent(inseq)))
	metrics.append("Match: "+ verifystructure(inseq,instruct))
	metrics.append("Average Structure Edit Distance: " + str(EditDistanceList[counter]))

	return metrics
def RetrieveFile(number_of_sequences,structure_and_sequence_list):
	global EditDistanceList,seqnumber
	seqnumber=number_of_sequences
	#print(len(structure_and_sequence_list))
	targetstructs= structure_and_sequence_list[0::2]
	counter=0
	structures= structure_and_sequence_list[0::2]
	print("Finding edit Dist...")
	EditDistanceList=Levenshtein_distance(structures)
	print("Generating Metrics...")
	for item in structure_and_sequence_list[1::2]:
		metriclist.append(GenMetrics(item,targetstructs[counter],counter))
		counter+=1
	#print(counter)

	#InverseSequenceOutputFile.close()
def GenMetricFile(gkm_svm_output):
	f=open("SequenceMetrics.txt", "w")
	iterable=0
	for item in metriclist:
		for item2 in item:
			#print(item2)
			f.write(item2)
			f.write('\n')
		f.write("normality: "+str(gkm_svm_output[iterable]))
		iterable+=1
		f.write('\n')
		f.write('\n')
	f.close()





