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

import os
import random
import shutil
import subprocess
import time
from multiprocessing import Pool
import SequenceMetricGeneration
import pandas as pd
import re
import itertools
import math
import PostGenerationMetrics

millis = int(round(time.time() * 1000))


# run a command and return the stdout from it
def stdout_from_command(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    return iter(p.stdout.readline, b'')


# RNAfold courtesy viennaRNA
def file_len(fname):
    with open(fname) as g:
        for i, l in enumerate(g):
            pass
    return i + 1


def get_moststable(sequence):  # RNALfold Stable Structure Finder
    a = []
    if len(sequence) < 1:
        return (0, 0)
    lines = stdout_from_command("echo %s | RNALfold --span=50" % sequence)
    for i in range(200):  # Get all the outputs. Probably not over 200. If it is, make this number bigger rinse repeat
        try:
            p = lines.next()
            a.append(p)
        except:
            pass
    stableseq = (a[-2])[0:-1]
    stableenergy = (float((a[-1])[2:-2]))
    return stableseq, stableenergy


def gen_diversity():
    with open('gkm_svm_output.txt') as j:
        line_count = 0
        for line in j:
            line_count += 1
    j.close()
    diversityscores = []
    f = open("gkm_svm_output.txt", "r")
    for z in range(0, line_count):
        temp = f.readline()
        if temp != "" and temp != "\n":
            if temp[1] == "1":
                sub = temp.split()
                diversityscores.append(sub[1])
    return diversityscores


# legend - numbering around the start codon

#  -6 -5 -4 -3 -2 -1 +1 +2 +3 +4 +5 
#  N  N  N  N  N  N  A  U  G  N  N

# # start codons to use
# start_codons = ["AUG", "CUG", "GUG", "UUG"]

# # stop codons
# stop_codons = ["UAA", "UGA", "UAG"] 

# start codons to use
start_codons = ["aug", "cug", "gug", "uug"]

# stop codons
stop_codons = ["uaa", "uga", "uag"]

# start codon contexts to use 
# start_context = ["GCCACC", "GCCGCC", "GUUACC", "GUUGCC", "GUUUCC"]

# -6 to -1 position start codon contexts from Noderer et al 2014 (weak to strong) 
start_context_pm6pm1 = ["UGAUAU", "AUCUUC", "UGCUUG", "GGCGCU", "AGAGUG", "AGCGCU", "UGUGGA", "UGCGUG", "AUCGCA"]

# +4 and +5 start codon contexts from Noderer (weak to strong)
# start_context_p4p5 = ["CC", "AU", "UC", "CA", "AC", "GC"]
start_context_p4p5 = ["cc", "au", "uc", "ca", "ac", "gc"]
# all distances in nucleotides (not codons) 

# distance between uORF start and structure - 8nt is mostly upstream, 14nt is ~100% upstream, and 32 nt is ~50/50 (Kozak PNAS 1990)
#  (above with uAUG context uuuAUGg and a stem loop of sequence GATCCGGGTTCTCCCGGATC (-14.2 kcal/mol) 
dist_uorf_strx = []

# distance between uORF stop and start codon
# - spacing between uORF stop and main ORF start (low -> high reinitiation: 11, 45, 79 nt between stop and downstream AUG; Kozak 1987 MCB)
dist_uorf_stop_main_start = [10, 20, 30, 40, 50]

# uORF length (not including start or stop codons) 
uorf_length_min = [0, 3, 12, 24, 39, 57, 99]


InverseReady = []
startcodonstartlocations=[]
insertrange=[]
stem_length = []
loop_length = []

def bindpattern(inpat, inseq):  # tool to bind the dot-bracket notation to its inherent pattern for RNAinverse
    return inpat + "\n" + inseq


def insert_with_delete(base, insert, index):  # tool to insert a string into another string at index with delete
    return base[0:index] + insert + base[index + len(insert):]


def insert(inputstring, index, insertedstr):  # tool to insert a string into another string at index without delete
    return inputstring[:index] + insertedstr + inputstring[index + 1:]


def generate_nucleotide_sequence(length,locs):  # writes a sequence of uORF length ulen and NCR length distusms intended for the nucleotide part of an RNAInverse call #TODO: Add negative context permutations

    retlist = []
    start_codons=['aug']
    for startcodon in start_codons:
        for location in locs:
            retlist.append(location*'N'+(startcodon + ((length-3-location) * 'N' )))
    return retlist



def takeinversefromstring(instring):  # simplified functionality for the inverse call from the input string

    uname = str(random.randint(0, 1000000))
    open(uname + ".txt", "w").close()
    file1 = open(uname + ".txt", "w")
    tmp = stdout_from_command("echo \'%s\' | RNAinverse -Fmp -f 0.5 -d2" % instring)

    lead = tmp.__next__().decode("utf-8")
    seq = tmp.__next__().split()[0].decode("utf-8")
    print(instring)
    print(lead.rstrip())
    print(seq + "\n")
    file1.write(instring + '\n' + lead.rstrip() + '\n' + seq + '\n\n')
    file1.close()
    return seq


def fullresource_inverse(adjoinedlist):
    p = Pool()  # start multi-core processing pool with all available resources
    result = p.map(takeinversefromstring,
                   adjoinedlist)  # use the pool toward processing the sequential list through the RNAinverse program
    p.close()  # stop the pool
    p.join()
    print("Done!")
    # print(result)


def length_based_generation(length, permutationlist):
    generatedseqs = []
    pairrdy = []

    # generate all the permutations of the nucleotide seq. Note: Constant offset vals due to feature (codons,
    # contexts) lengths.


    generatedseqs.append(generate_nucleotide_sequence(length,startcodonstartlocations))
    print(generatedseqs)
    orderednts = [item for sublist in generatedseqs for item in sublist]
    print(orderednts)
    print(len(permutationlist), "structural permutations")
    print(len(orderednts), "nucleotide permutations")

    start = time.time()
    # 1 nucleotide sequence to 1 structure
    for nucleotidesequence in orderednts:
        for permutation in permutationlist:
            pair = [nucleotidesequence, permutation]
            pairrdy.append(pair)

    for item in pairrdy:
        InverseReady.append(bindpattern(item[1], item[0]))

    end = time.time()


    print(len(InverseReady), "inverse ready pairs")


# assembles the individually generated txt files into a unified text file
def assemble():
    global structs_and_seqs
    namelist = os.listdir(project_path + 'AssemblySeq/')
    sequenceoutputlist = []
    for name in namelist:
        with open(name) as f:
            with open('InverseSequenceOutput.txt', 'a') as f1:
                for line in f:
                    f1.write(line)
                    sequenceoutputlist.append(line)
    try:
        f1.close()
    except:
        pass

    shutil.copyfile(project_path + 'AssemblySeq/InverseSequenceOutput.txt',
                    project_path + 'InverseSequenceOutput.txt')
    os.chdir(project_path)
    shutil.rmtree('AssemblySeq')
    # print(sequenceoutputlist)
    structs = [x[:-1] for x in sequenceoutputlist[0::5]]
    # print(structs)
    seqs = [x[:-1] for x in sequenceoutputlist[3::5]]
    # print(seqs)
    structs_and_seqs = [x for y in zip(structs, seqs) for x in y]


def fold(sequence):  # From Original Program
    lines = stdout_from_command("echo %s | RNAfold --noPS --MEA -p" % sequence)

    # first line is just the seq - skip
    lines.__next__()
    secstruct = lines.__next__().decode("utf-8").split()[0]
    return secstruct


def generate_summary_data(sequence_count):  # Previously in PostGenerationMetrics.py. Returns library summary data
    v = open("SequenceMetrics.txt", "r")
    start_position_1 = 11
    start_position_2 = 2
    start_position_3 = 3
    listofseqs = []
    listofnormvals = []
    triplelist = []
    structlib = []
    verificationlist = []
    a = 0
    iter = 0
    for x in range(12 * sequence_count):
        try:
            temp = v.readline()

            if (temp[7:10] == "Ver" or temp[7:10] == "Dis"):
                verificationlist.append(temp[7:-1])
            elif (temp[0] == "." or temp[0] == "("):
                structlib.append(temp[0:-1])

            elif (start_position_1 % 12 == 0):  # the sequence
                triplelist.append(temp[0:-1])

            elif (start_position_2 % 12 == 0):  # the normality
                triplelist.append(float(temp[11:-1]))
            elif (start_position_3 % 12 == 0):  # the normality
                triplelist.append(float(temp[33:-1]))

            start_position_1 += 1
            start_position_2 += 1
            start_position_3 += 1
            iter += 1
        except:
            pass
    verified_count = 0
    for item in verificationlist:
        if (item[0] == 'V'):
            verified_count += 1
    print('Percent Verified: ' + str(verified_count / sequence_count * 100))

    return verificationlist, structlib, triplelist[0::3], triplelist[1::3], triplelist[2::3]


# specifies the path of where you want to run this script. Note that if you are running
# the R script or Metrics, those files have to be in the same directory
project_path = '/Users/Peter/PycharmProjects/RNALibrary/'


def generate(sequence_length, permutationlist):  # the runs the generation of the sequences
    # all of this just clears the files for a re-run

    os.chdir(project_path)
    open('InverseSequenceOutput.txt', 'w').close()
    os.chdir(project_path)

    try:
        os.mkdir(project_path + 'AssemblySeq/')
    except:
        shutil.rmtree('AssemblySeq')
        os.mkdir(project_path + 'AssemblySeq/')

    os.chdir(project_path + 'AssemblySeq/')

    # generation and inverse of the sequences
    length_based_generation(sequence_length, permutationlist)
    # print(InverseReady[0:sequence_count])
    fullresource_inverse(InverseReady)
    assemble()

    subprocess.call("Rscript gkm_svm_generator.R", shell=True)
    gkm_svm_scores = gen_diversity()

    SequenceMetricGeneration.RetrieveFile(len(InverseReady), structs_and_seqs)
    SequenceMetricGeneration.GenMetricFile(gkm_svm_scores)


structs_and_seqs = []

data = pd.read_csv(r'/Users/Peter/Downloads/snv_phenotype_log_diff.csv')
df = pd.DataFrame(data, columns=['utr', 'info4']).to_numpy()
nonvariants = []
for sequencevariantpair in df:
    if sequencevariantpair[1] == "normal":
        nonvariants.append(sequencevariantpair[0].replace("T", "U"))


def generate_pin(stemlen, looplen, insertnumber):  # generate a hairpin structure in dot-bracket notation

    base = ("(" * stemlen + "." * looplen + ")" * stemlen)

    for item in range(insertnumber):
        decider = random.randint(0, 1)
        if (decider == 1):
            randindex = random.randint(1, stemlen)
        else:
            randindex = random.randint(stemlen + looplen + 1, stemlen * 2 + looplen - 1)
        base = base[:randindex] + "." + base[randindex:]
    return base


def get_moststable(sequence,seqnumber):  # RNALfold Stable Structure Finder
    a = []
    if len(sequence) < 1:
        return (0, 0)
    lines = stdout_from_command("echo %s | RNALfold --span=50" % sequence)
    for i in range(
            seqnumber * 2):  # Get all the outputs. Probably not over 200. If it is, make this number bigger rinse repeat
        # print(a)
        try:
            a.append(lines.__next__())
        # print(a)
        except:
            pass
    stableseq = (a[-2])[0:-1]
    stableenergy = (float((a[-1])[2:-2]))
    return stableseq, stableenergy


def generate_structure_permutations(nucleotide_sequence):  # generate all structural permutations within the
    # possibility of the given parameters stem length and loop length. Throws away sequences if impossible
    permutations = []
    advancedpermutations = []
    blankpattern = len(nucleotide_sequence) * "."




    random.shuffle(stem_length)
    random.shuffle(loop_length)
    random.shuffle(insertrange)




    for slen in stem_length:
        for llen in loop_length:
            for dus in dist_uorf_strx:
                for insertnum in insertrange:
                    permutation = insert_with_delete(blankpattern, generate_pin(slen, llen, insertnum), dus)
                    if (len(permutation) <= len(
                            blankpattern)):  # throws out structures which are larger than possible within the configuration
                        permutations.append(permutation)
                        # print(permutation)

    return permutations





if __name__ == '__main__':
    #global startcodonstartlocations

    #Depends on 7 variables: sequence_length, structuralcount, locs (product of these the first and the length of the second = seq #),
    #dist_uorf_strx, insertrange, loop_length, and stem_length (change the features of the structure)

    sequence_length = 30

    structuralcount = 2
    startcodonstartlocations=[2,9,17,24,27]

    dist_uorf_strx = [8] #where the first start codon is
    insertrange = [0 ,2 ,3] #how many dots within the hairpin other than the loop
    stem_length = [3, 4, 5]
    loop_length = [3, 4, 5]


    permutationlist = []
    start = time.time()


    permutationlist = random.sample(generate_structure_permutations("n" * sequence_length),structuralcount)

    generate(sequence_length, permutationlist)

