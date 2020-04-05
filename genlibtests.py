#!/usr/bin/env python

# program to generate library of RNA sequences with defined combinations of translation regulatory elements

# elements to include / variables: 
# - out of frame uAUG
# - uORF (uAUG with stop codon)
# - spacing between uORF stop and main ORF start (low -> high reinitiation: 11, 45, 79 nt between stop and downstream AUG; Kozak 1987 MCB)
# - structure folding free energy/hairpin length
# - distance between structure(s) and uAUG of uORF
# - main AUG kozak
# - uAUG kozak
# - codon preceding uORF stop codon (Grant & Hinnebusch 1994)
# - uORF length (up to 39 nt is ok, 50% efficiency of reinitiation after 57 nt; Hwang & Su JGV 1998 or Child & Geballe 1999)

# assumptions
# - according to a meeting with Twist, they can synthesize most anything in the oligo pool format (can send them the sequences and they can check)



# depends: RNAinverse (ViennaRNA; externally installed and in PATH) 

# TODO:

# - convert U to T in output
# - add some sort of handle/etc for cloning, filter sequences if necessary 
# - sumbout 5' TOP 
# - distance btw 5' to uAUG - can enhance initiation when 5' leader is > 35 in yeast (Slusher et al 1991) 
#testcommit
import sys, os, argparse, subprocess, re, threading, numpy, time

millis=int(round(time.time()*1000))
# run a command and return the stdout from it 
def stdout_from_command(command):
    p = subprocess.Popen(command, stdout = subprocess.PIPE, shell = True)
    return iter(p.stdout.readline, b'')


# gets the folding energy of a sequence. returns a tuple containing energies (centroid_energy, structure)
# RNAfold courtesy viennaRNA
def get_rna_structure_energy(sequence):
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

    return (nrg, strx)
                                                                                                                          
# legend - numbering around the start codon

#  -6 -5 -4 -3 -2 -1 +1 +2 +3 +4 +5 
#  N  N  N  N  N  N  A  U  G  N  N

# start codons to use
start_codons = ["AUG, CUG, GUG, UUG"]

# stop codons
stop_codons = ["TAA", "TGA", "TAG"] 

# start codon contexts to use 
#start_context = ["GCCACC", "GCCGCC", "GUUACC", "GUUGCC", "GUUUCC"]

# -6 to -1 position start codon contexts from Noderer et al 2014 (weak to strong) 
start_context_pm6pm1 = ["UGAUAU", "AUCUUC", "UGCUUG", "GGCGCU", "AGAGUG", "AGCGCU", "UGUGGA", "UGCGUG", "AUCGCA"]

# +4 and +5 start codon contexts from Noderer (weak to strong)
start_context_p4p5 = ["CC", "AU", "UC", "CA", "AC", "GC"]

# all distances in nucleotides (not codons) 

# distance between uORF start and structure - 8nt is mostly upstream, 14nt is ~100% upstream, and 32 nt is ~50/50 (Kozak PNAS 1990)
#  (above with uAUG context uuuAUGg and a stem loop of sequence GATCCGGGTTCTCCCGGATC (-14.2 kcal/mol) 
dist_uorf_strx = [8, 12, 14, 16, 20, 24, 32]

# distance between uORF stop and start codon
# - spacing between uORF stop and main ORF start (low -> high reinitiation: 11, 45, 79 nt between stop and downstream AUG; Kozak 1987 MCB)
dist_uorf_stop_main_start = [10, 20, 30, 40, 50]

# uORF length (not including start or stop codons) 
uorf_length = [0, 3, 12, 24, 39, 57, 99]

stem_length = [4]

loop_length= [6]


max_lib_size = len(start_codons) * len(stop_codons) * len(start_context_pm6pm1) * len(start_context_p4p5) * \
               len(dist_uorf_strx) * len(dist_uorf_stop_main_start) * len(uorf_length)

rnainverse_input_long = "...(((((....)))))....(((.(((....))).))).......(((((....)))))....(((.(((....))).))).......(((((....)))))....(((.(((....))).))).......(((((....)))))....(((.(((....))).))).......(((((....)))))....(((.(((....))).)))....\naugNNNNNugauNNNNNaugaNNNuNNNguuaNNNuNNNNNNNaugNNNNNugauNNNNNaugaNNNuNNNguuaNNNuNNNNNNNaugNNNNNugauNNNNNaugaNNNuNNNguuaNNNuNNNNNNNaugNNNNNugauNNNNNaugaNNNuNNNguuaNNNuNNNNNNNaugNNNNNugauNNNNNaugaNNNuNNNguuaNNNuNNNNNNN"

rnainverse_input_short = "...(((((....)))))....(((.(((....))).)))...\naugNNNNNugauNNNNNaugaNNNuNNNguuaNNNuNNNNNN"


a="."
b="("
c=")"
masterlist=[]
def bindpattern(inseq,inpat): #tool to bind the dot-bracket notation to its inherent pattern
    return inseq+"\n"+inpat

def takeinverse(instring): #simplified functionality for the inverse call from the input string
    tmp = stdout_from_command("echo \'%s\' | RNAinverse -Fmp -f 0.5 -d2" % instring)
    tmp.next()
    seq = tmp.next().split()[0]
    print(seq)
    return seq
def insert_with_delete(base,insert,index):#tool to insert a string into another string at index with delete
    return base[0:index]+insert+base[index+len(insert):]

def insert(inputstring, index, insertedstr): #tool to insert a string into another string at index without delete
    return inputstring[:index]+insertedstr+inputstring[index+1:]
    
def genpin(pinlen,looplen): #generate a hairpin structure in dot-bracket notation
    return(b*pinlen+a*looplen+c*pinlen)

def generate_stop_start_pattern(ulen,distusms): #writes a sequence of uORF length ulen and NCR length distusms intended for the latter part of an RNAInverse call #TODO: Add negative context permutations
    return "aug"+ulen*"N"+"uag"+(distusms)*"N"+"aug"

def generate(): #still developing function which will serve as the primary generation. Will iterate through the input parameters, generating and inverse folding all permutation
    for lengths in uorf_length:
        for distances in dist_uorf_stop_main_start:
            ntpattern=generate_stop_start_pattern(lengths,distances) #generates the defined nucleotides and undefined regions of the sequence
            basepat=(lengths+6)*a+(distances+3)*a #generates a "blank" pattern for the length of the structural sequence
            for loops in loop_length:
                for stems in stem_length:
                    pin=genpin(loops,stems) #generate the hairpin structure 
                    for startx in dist_uorf_strx:
                        patwithstruct=insert_with_delete(basepat,pin,startx) #insert with delete the hairpin structure into the structural sequence
                        boundpat=bindpattern(ntpattern,patwithstruct) #concatenate the nucleotide pattern to the structural pattern
                        masterlist.append(boundpat)

     
    
    #inverse=takeinverse(boundpat) #uses RNAinverse to take the inverse of the given nucleotide pattern and structural sequence
    return #boundpat
print("total sequences "+str(len(uorf_length)*len(dist_uorf_stop_main_start)*len(loop_length)*len(stem_length)*len(dist_uorf_strx)))
generate()
#print(masterlist)
masterlist1=masterlist[0:len(masterlist)/2]
masterlist2=masterlist[(len(masterlist)/2):]
masterinverse=[]


def takeinversefromlist(instring): #simplified functionality for the inverse call from the input string
    
    tmp = stdout_from_command("echo \'%s\' | RNAinverse -Fmp -f 0.5 -d2" % instring)
    tmp.next()
    seq = tmp.next().split()[0]
    print(seq)
    masterinverse.append(seq)
    return



#numpy.savetxt("output.csv", a, delimiter=",")
#a=numpy.append(a,[3,5], axis=1)
print(a)
def inverse_with_multithreading(adjoinedlist):

    masterlist1=adjoinedlist[0:len(adjoinedlist)/2]
    masterlist2=adjoinedlist[(len(adjoinedlist)/2):]
    def thread1():
        timeseries_thread1 = numpy.array([[0, 0]])
        k=0
        for thread1iter in masterlist1:
            thread1_basetime=int(round(time.time()*1000))
            takeinversefromlist(thread1iter)
            #print("thread 1: "+str(k)+" of "+str((len(masterlist1))))
            thread1_computetime=(int(round(time.time()*1000)-thread1_basetime))
            #print("thread 1 compute time: "+str(thread1_computetime))
            thread1append= numpy.array([[k,thread1_computetime]])
            timeseries_thread1=numpy.append(timeseries_thread1,thread1append,axis=0)
            k+=1

        numpy.savetxt("thread1.csv", timeseries_thread1, delimiter=",")
    def thread2():
        timeseries_thread2 = numpy.array([[0, 0]])
        v=(len(adjoinedlist)/2)
        for thread2iter in masterlist2:
            thread2_basetime=int(round(time.time()*1000))
            takeinversefromlist(thread2iter)
            #print("thread 2: "+str(v)+" of "+str(len(masterlist)))
            thread2_computetime=(int(round(time.time()*1000)-thread2_basetime))
            #print("thread 2 compute time: "+str(thread2_computetime))
            thread2append= numpy.array([[v,thread2_computetime]])
            timeseries_thread2=numpy.append(timeseries_thread2,thread2append,axis=0)
            v+=1
        numpy.savetxt("thread2.csv", timeseries_thread2, delimiter=",")
    
    t1 = threading.Thread(target=thread1, args=())
    t2 = threading.Thread(target=thread2, args=())
    t1.start() 
    # starting thread 2 
    t2.start()
    # wait until thread 2 is completely executed
    t1.join() 
    # wait until thread 2 is completely executed 
    t2.join() 
    
  
    # starting thread 1 
    
  
    # both threads completely executed 
    print("Done!") 
inverse_with_multithreading(masterlist)



#dist_uorf_strx = [8, 12, 14, 16, 20, 24, 32]

# distance between uORF stop and start codon
# - spacing between uORF stop and main ORF start (low -> high reinitiation: 11, 45, 79 nt between stop and downstream AUG; Kozak 1987 MCB)
#dist_uorf_stop_main_start = [10, 20, 30, 40, 50]

# uORF length (not including start or stop codons) 
#uorf_length = [0, 3, 12, 24, 39, 57, 99]







