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
dist_uorf_strx = list(range(0, 35))

# distance between uORF stop and start codon
# - spacing between uORF stop and main ORF start (low -> high reinitiation: 11, 45, 79 nt between stop and downstream AUG; Kozak 1987 MCB)
dist_uorf_stop_main_start = [10, 20, 30, 40, 50]

# uORF length (not including start or stop codons) 
uorf_length_min = [0, 3, 12, 24, 39, 57, 99]

stem_length = [3, 4, 5]

loop_length = [3, 4, 5]
InverseReady = []


def bindpattern(inpat, inseq):  # tool to bind the dot-bracket notation to its inherent pattern for RNAinverse
    return inpat + "\n" + inseq


def insert_with_delete(base, insert, index):  # tool to insert a string into another string at index with delete
    return base[0:index] + insert + base[index + len(insert):]


def insert(inputstring, index, insertedstr):  # tool to insert a string into another string at index without delete
    return inputstring[:index] + insertedstr + inputstring[index + 1:]


def generate_pin(stemlen, looplen):  # generate a hairpin structure in dot-bracket notation
    return ("(" * stemlen + "." * looplen + ")" * stemlen)


def generate_nucleotide_sequence(ulen,
                                 distusms):  # writes a sequence of uORF length ulen and NCR length distusms intended for the nucleotide part of an RNAInverse call #TODO: Add negative context permutations
    retlist = []
    for elem in start_codons:
        for elem2 in start_codons:
            for elem3 in stop_codons:
                for contexts in start_context_p4p5:
                    pt1l3 = elem
                    pt2l3 = elem3
                    pt3l8 = contexts + 3 * "N" + elem2
                    retlist.append((pt1l3 + ulen * 'N' + pt2l3 + distusms * 'N' + pt3l8))
    return retlist


def generate_structure_permutations(nucleotide_sequence):  # generate all structural permutations within the
    # possibility of the given parameters stem length and loop length. Throws away sequences if impossible
    permutations = []
    blankpattern = len(nucleotide_sequence) * "."
    for slen in stem_length:
        for llen in loop_length:
            for dus in dist_uorf_strx:
                permutation = insert_with_delete(blankpattern, generate_pin(slen, llen), dus)
                if (len(permutation) <= len(
                        blankpattern)):  # throws out structures which are larger than possible within the configuration
                    permutations.append(permutation)
    return permutations


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
    #print(result)


def length_based_generation(length):
    generatedseqs = []
    pairrdy = []

    # generate all the permutations of the nucleotide seq. Note: Constant offset vals due to feature (codons,
    # contexts) lengths.
    for i in range(0, length - 14 + 1):
        generatedseqs.append(generate_nucleotide_sequence(i, length - 14 - i))

    orderednts = [item for sublist in generatedseqs for item in sublist]

    for nucleotidesequence in orderednts:
        for structuralsequence in generate_structure_permutations(nucleotidesequence):
            pair = [nucleotidesequence, structuralsequence]
            pairrdy.append(pair)

    for item in pairrdy:
        InverseReady.append(bindpattern(item[1], item[0]))


    # for library diversity given limited computational power. if all sequences are generated then this is not needed
    random.shuffle(InverseReady)
    print(len(InverseReady))


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
    verified_count=0
    for item in verificationlist:
        if(item[0]=='V'):
            verified_count+=1
    print('Percent Verified: ' +str(verified_count/sequence_count*100))

    return verificationlist, structlib, triplelist[0::3], triplelist[1::3], triplelist[2::3]


# specifies the path of where you want to run this script. Note that if you are running
# the R script or Metrics, those files have to be in the same directory
project_path = '/Users/Peter/PycharmProjects/RNALibrary/'


def generate(sequence_count, sequence_length):  # the runs the generation of the sequences
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
    length_based_generation(sequence_length)
    fullresource_inverse(InverseReady[0:sequence_count])
    assemble()

    subprocess.call("Rscript gkm_svm_generator.R", shell=True)
    gkm_svm_scores = gen_diversity()

    SequenceMetricGeneration.RetrieveFile(sequence_count, structs_and_seqs)
    SequenceMetricGeneration.GenMetricFile(gkm_svm_scores)


structs_and_seqs = []
if __name__ == '__main__':
    sequence_count = 30
    sequence_length = 45
    generate(sequence_count, sequence_length)

    verificationlist, structlib, seqs, edit, diversity = generate_summary_data(sequence_count)


