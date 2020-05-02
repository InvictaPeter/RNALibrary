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

import sys, os, argparse, subprocess, re, multiprocessing, numpy, time
from multiprocessing import Pool
from subprocess import Popen
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
start_codons = ["AUG", "CUG", "GUG", "UUG"]

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
uorf_length_min = [0, 3, 12, 24, 39, 57, 99]

stem_length = [4]

loop_length= [3]



def bindpattern(inseq,inpat): #tool to bind the dot-bracket notation to its inherent pattern
    return inseq+"\n"+inpat

def insert_with_delete(base,insert,index):#tool to insert a string into another string at index with delete
    return base[0:index]+insert+base[index+len(insert):]

def insert(inputstring, index, insertedstr): #tool to insert a string into another string at index without delete
    return inputstring[:index]+insertedstr+inputstring[index+1:]
    
def generate_pin(stemlen,looplen): #generate a hairpin structure in dot-bracket notation
    return("("*stemlen+"."*looplen+")"*stemlen)


def generate_nucleotide_sequence(ulen,distusms): #writes a sequence of uORF length ulen and NCR length distusms intended for the nucleotide part of an RNAInverse call #TODO: Add negative context permutations
    retlist=[]
    for elem in start_codons:
        for elem2 in start_codons:
            for elem3 in stop_codons:
                for contexts in start_context_p4p5:
                    retlist.append((elem+ulen*"N"+elem3+(distusms-5)*"N"+contexts+3*"N"+elem2))
    return retlist

def generate_nucleotide_sequence_RD(ulen,distusms): #writes a sequence of uORF length ulen and NCR length distusms intended for the nucleotide part of an RNAInverse call #TODO: Add negative context permutations
    retlist=[]
    for elem in start_codons:
        for elem2 in start_codons:
            for elem3 in stop_codons:
                for contexts in start_context_p4p5:
                    pt1l3=elem
                    pt2l3=elem3
                    pt3l8=contexts+3*"N"+elem2
                    #retlist.append((elem+(ulen-6)*"N"+elem3+(distusms-5-6)*"N"+contexts+3*"N"+elem2))
                    retlist.append((pt1l3+ulen*'N'+pt2l3+distusms*'N'+pt3l8))
    return retlist

def generate_structure_permutations(nucleotide_sequence): #generate all structural permutations within the possibility of the given parameters stem length and loop length. Throws away sequences if impossible
    permutations=[]
    blankpattern=len(nucleotide_sequence)*"."
    for slen in stem_length:
        for llen in loop_length:
            for dus in dist_uorf_strx:
                permutation=insert_with_delete(blankpattern,generate_pin(slen,llen),dus)
                if(len(permutation)<=len(blankpattern)): #throws out structures which are larger than possible within the configuration
                    permutations.append(permutation)
    return permutations

def generate_nucleotide_permutations(): #generates the nucleotide (nonstructural) permutations within given parameters.
    nucleotide_permutations=[]
    for distances in dist_uorf_stop_main_start:
        for lengths in uorf_length:
            for sequences in generate_nucleotide_sequence(lengths,distances):
                nucleotide_permutations.append(sequences)
    nucleotide_permutations.sort(key=len)
    return nucleotide_permutations

def takeinversefromstring(instring): #simplified functionality for the inverse call from the input string
    file1 = open("subbatch1.txt","a") 
    global bruh
    tmp = stdout_from_command("echo \'%s\' | RNAinverse -Fmp -f 0.5 -d2" % instring)
    
    lead=tmp.next()
    seq = tmp.next().split()[0]
    print(instring)
    print(lead.rstrip())
    print(seq+"\n")
    file1.write(instring+'\n'+lead.rstrip()+'\n'+seq+'\n\n')
    return seq
def inverse_with_multithreading(adjoinedlist):
    p=Pool() #start multi-core processing pool with all available resources
    result=p.map(takeinversefromstring,adjoinedlist) #use the pool toward processing the sequential list through the RNAinverse program
    p.close() #stop the pool
    p.join()
    print("Done!") 
    print(result)

#def main():
def LengthGen(length):
    a=[]
    g=[]
    orderednts=[]
    pairrdy=[]
    InverseReady=[]
    for i in range(0,length-14+1): #generate all the permutations of the nucleotide seq. NOTE TO SELF: DONT PLAY WITH THE CONSTANTS HERE. THEY WORK.
        g.append(generate_nucleotide_sequence_RD(i,length-14-i))
    # print((g[0])[0])
    # print((g[1])[0])
    # print((g[2])[0])
    for e in range(len((g[0])[0])):
        for i in range(len(g)):
            orderednts.append((g[i])[e])
    for nucleotidesequence in orderednts:
        for structuralsequence in generate_structure_permutations(nucleotidesequence):
            pair=[nucleotidesequence,structuralsequence]
            pairrdy.append(pair)
    #print(pairrdy)
    for item in pairrdy:
        InverseReady.append(bindpattern(item[0],item[1]))
    print(len(InverseReady))
    # for i in generate_structure_permutations((g[0])[0]):
    #     print(len(i))
    # for i in range(len(g)):
    #     print((g[i])[0])

    
if __name__ == '__main__':
    LengthGen(150)
    #open('subbatch1.txt', 'w').close()
    #main()



