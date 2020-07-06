#!/usr/bin/env python
"""
Created on Apr 10, 2015
Modified for RNA-seq on Apr 8, 2018
Modified again in June 2018
@first_author: Chen Yang
@secondary_author: Lucus Pocus

This script generates simulated Oxford Nanopore 1D RNA-seq reads from one reference transcript sequence.

"""


from __future__ import print_function
from __future__ import with_statement
import sys
import getopt
import random
import re
from time import strftime
try:
    from six.moves import xrange
except ImportError:
    pass
try:
    import numpy as np
except ImportError:
    sys.exit("""You need numpy!
                install it from http://www.numpy.org/""")
import mixed_model as mm

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("""You need BioPython!
                install it from ...""")
    
PYTHON_VERSION = sys.version_info
VERSION = "1.0.0"
PRORAM = "RNanoSim"

BASES = ['A', 'T', 'C', 'G']


# Usage information
def usage():
    usage_message = "./transimulator.py [command] <options>\n" \
                    "<options>: \n" \
                    "-h : print usage message\n" \
                    "-r : reference transcripts sequences in fasta file, specify path and file name, REQUIRED\n" \
                    "--select : file containing the id of transcripts to be simulated (should match those in reference fasta file)\n" \
                    "-c : The prefix of training set profiles, same as the output prefix in read_analysis.py, default = training\n" \
                    "-o : The prefix of output file, default = 'simulated'\n" \
                    "-n : Number of generated reads, default = 50 reads\n" \
                    "--KmerBias: prohibits homopolymers with length >= n bases in output reads, default = 6\n" \
                    "--seed: manually seeds the pseudo-random number generator, default = None\n"

    sys.stderr.write(usage_message)






###############################################################################
# FUNCTION read_ecdf()
#    
###############################################################################
    
def read_ecdf(profile):
    # We need to count the number of zeros. If it's over 10 zeros, l_len/l_ratio need to be changed to higher.
    # Because it's almost impossible that the ratio is much lower than the lowest heuristic value.
    header = profile.readline()
    header_info = header.strip().split()
    ecdf_dict = {}
    lanes = len(header_info[1:])

    for i in header_info[1:]:
        boundaries = i.split('-')
        ecdf_dict[(int(boundaries[0])), int(boundaries[1])] = {}

    ecdf_key = sorted(ecdf_dict.keys())
    l_prob = [0.0] * lanes
    l_ratio = [0.0] * lanes

    for line in profile:
        new = line.strip().split('\t')
        ratio = [float(x) for x in new[0].split('-')]
        prob = [float(x) for x in new[1:]]
        for i in xrange(lanes):
            if prob[i] == l_prob[i]:
                continue
            else:
                if l_prob[i] != 0:
                    ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] = (l_ratio[i], ratio[1])
                else:
                    ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] \
                        = (max(l_ratio[i], ratio[1] - 10 * (ratio[1] - ratio[0])), ratio[1])
                l_ratio[i] = ratio[1]
                l_prob[i] = prob[i]

    for i in xrange(0, len(ecdf_key)):
        last_key = sorted(ecdf_dict[ecdf_key[i]].keys())[-1]
        last_value = ecdf_dict[ecdf_key[i]][last_key]
        ecdf_dict[ecdf_key[i]][last_key] = (last_value[0], ratio[1])

    return ecdf_dict

# Original: stochastically defines read length based on empirical distribution
# Modification: constant length
def get_length(len_dict, num, max_l, min_l):
    length_list = []
    for i in xrange(num):  
        middle_ref = max_l
        key = tuple(len_dict.keys())[0]
        while middle_ref <= min_l or middle_ref > max_l:
            p = random.random()
            for k_p, v_p in len_dict[key].items():
                if k_p[0] <= p < k_p[1]:
                    middle_ref = int(round((p - k_p[0])/(k_p[1] - k_p[0]) * (v_p[1] - v_p[0]) + v_p[0]))
                    break
        length_list.append(middle_ref)

    return length_list


def trans_read_profile(model_prefix):
    global out_reads
    global match_ht_list, align_ratio, ht_dict, error_par
    global trans_error_pr, match_markov_model

    # Read model profile for match, mismatch, insertion and deletions
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read error profile\n")
    sys.stdout.flush()
    error_par = {}
    model_profile = model_prefix + "_model_profile"
    with open(model_profile, 'r') as mod_profile:
        mod_profile.readline()
        for line in mod_profile:
            new_line = line.strip().split("\t")
            if "mismatch" in line:
                error_par["mis"] = [float(x) for x in new_line[1:]]
            elif "insertion" in line:
                error_par["ins"] = [float(x) for x in new_line[1:]]
            else:
                error_par["del"] = [float(x) for x in new_line[1:]]

    trans_error_pr = {}
    with open(model_prefix + "_error_markov_model", "r") as error_markov:
        error_markov.readline()
        for line in error_markov:
            info = line.strip().split()
            k = info[0]
            trans_error_pr[k] = {}
            trans_error_pr[k][(0, float(info[1]))] = "mis"
            trans_error_pr[k][(float(info[1]), float(info[1]) + float(info[2]))] = "ins"
            trans_error_pr[k][(1 - float(info[3]), 1)] = "del"

    with open(model_prefix + "_first_match.hist", 'r') as fm_profile:
        match_ht_list = read_ecdf(fm_profile)

    with open(model_prefix + "_match_markov_model", 'r') as mm_profile:
        match_markov_model = read_ecdf(mm_profile)

    # Read profile of aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read ECDF of aligned reads\n")
    sys.stdout.flush()

    # Read align ratio profile
    # with open(model_prefix + "_align_ratio", 'r') as a_profile:
    #    align_ratio = read_ecdf(a_profile)

    # Read head/unaligned region ratio
    with open(model_prefix + "_ht_ratio", 'r') as ht_profile:
        ht_dict = read_ecdf(ht_profile)
 
    #---> currently length is determined by reference transcript length.


def collapse_homo(seq, k):
    read = re.sub("A" * k + "+", "A" * (k - 1), seq)
    read = re.sub("C" * k + "+", "C" * (k - 1), read)
    read = re.sub("T" * k + "+", "T" * (k - 1), read)
    read = re.sub("G" * k + "+", "G" * (k - 1), read)

    return read


# Taken from https://github.com/lh3/readfq
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

#ref -> seq
#loop over ref sequences
#dna-type? remove?
#genome_len?
#def trans_simulation(ref, out, dna_type, per, kmer_bias, max_l, min_l):
def trans_simulation(seq_record, out, kmer_bias, number, out_reads, out_info, out_error):    
    global match_ht_list, ht_dict, match_markov_model
    global trans_error_pr, error_par
    
    seq_name = seq_record.id
    seq_len = len(seq_record.seq)
    
    i = 0
    while i < number:
        #ref = get_length(aligned_dict, 1, max_l, min_l)[0]
        ref = seq_len
        #simulate error profile (error_dict) and middle
        #middle: length of the "middle part" of the read that can be aligned on reference
        middle, middle_ref, error_dict = error_list(ref, match_markov_model, match_ht_list, error_par,
                                                    trans_error_pr)        
        new_read = seq_record.seq
        new_read_name = seq_name + "_" + str(i+1)
        
        # Mutate read
        new_read = case_convert(new_read)
        read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias)

        if kmer_bias:
            read_mutated = collapse_homo(read_mutated, kmer_bias)

        out_reads.write(">" + new_read_name + '\n')
        out_reads.write(read_mutated + '\n')
        
        out_info.write(new_read_name+"\t"+str(len(read_mutated))+"\t"+seq_name+"\t"+str(seq_len)+"\n")

        i += 1
        
    ht_dict.clear()


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq

#-> Needed??
def extract_read(dna_type, length):
    global seq_dict, seq_len, genome_len

    if length > max(seq_len.values()):
        length = max(seq_len.values())
    #*****
    # Extract the aligned region from reference
    if dna_type == "circular":
        ref_pos = random.randint(0, genome_len)
        chromosome = list(seq_dict.keys())[0]
        new_read_name = chromosome + "_" + str(ref_pos)
        if length + ref_pos <= genome_len:
            new_read = seq_dict[chromosome][ref_pos: ref_pos + length]
        else:
            new_read = seq_dict[chromosome][ref_pos:]
            new_read = new_read + seq_dict[chromosome][0: length - genome_len + ref_pos]
    else:
        # Generate a random number within the size of the genome. Suppose chromosomes are connected
        # tail to head one by one in the order of the dictionary. If the start position fits in one
        # chromosome, but the end position does not, then restart generating random number.
        while True:
            new_read = ""
            #ref_pos = random.randint(0, genome_len)
            #*****We want full-length transcripts
            ref_pos = 0
            for key in seq_len.keys():
                if ref_pos + length <= seq_len[key]:
                    new_read = seq_dict[key][ref_pos: ref_pos + length]
                    new_read_name = key + "_" + str(ref_pos)
                    break
                elif ref_pos < seq_len[key]:
                    break
                else:
                    ref_pos -= seq_len[key]
            if new_read != "":
                break
    return new_read, new_read_name


def unaligned_error_list(length, error_p):
    e_dict = {}
    error_rate = {(0, 0.4): "match", (0.4, 0.7): "mis", (0.7, 0.85): "ins", (0.85, 1): "del"}
    pos = 0
    last_is_ins = False
    while pos < length:
        p = random.random()
        for k_error in error_rate.keys():
            if k_error[0] <= p < k_error[1]:
                error_type = error_rate[k_error]
                break

        if error_type == "match":
            step = 1

        elif error_type == "mis":
            step = mm.pois_geom(error_p["mis"][0], error_p["mis"][2], error_p["mis"][3])
            e_dict[pos] = ["mis", step]

        elif error_type == "ins":
            step = mm.wei_geom(error_p["ins"][0], error_p["ins"][1], error_p["ins"][2], error_p["ins"][3])
            if last_is_ins:
                e_dict[pos + 0.1][1] += step
            else:
                e_dict[pos + 0.1] = ["ins", step]
                last_is_ins = True

        else:
            step = mm.wei_geom(error_p["del"][0], error_p["del"][1], error_p["del"][2], error_p["del"][3])
            e_dict[pos] = ["del", step]

        if error_type != "ins":
            pos += step
            last_is_ins = False

        if pos > length:
            length = pos

    return length, e_dict

###############################################################################
# FUNCTION error_list    
# INPUTs:
# m_ref: original read length (in RNanoSim -> reference transcript's length)   
# OUTPUTs:
# l_new:
# middle_ref:
# e_dict: dictionary read position -> [error type, step]   
###############################################################################   
def error_list(m_ref, m_model, m_ht_list, error_p, trans_p):
    # l_old is the original length, and l_new is used to control the new length after introducing errors
    l_new = m_ref
    pos = 0
    e_dict = {}
    middle_ref = m_ref
    prev_error = "start"

    # The first match come from m_ht_list
    p = random.random()
    k1 = list(m_ht_list.keys())[0]
    for k2, v2 in m_ht_list[k1].items():
        if k2[0] < p <= k2[1]:
            prev_match = int(np.floor((p - k2[0])/(k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
            if prev_match < 2:
                prev_match = 2
    pos += prev_match

    # Select an error, then the step size, and then a match and so on so forth.
    while pos < middle_ref:
        # pick the error based on Markov chain
        p = random.random()
        for k in trans_p[prev_error].keys():
            if k[0] <= p < k[1]:
                error = trans_p[prev_error][k]
                break

        if error == "mis":
            step = mm.pois_geom(error_p["mis"][0], error_p["mis"][2], error_p["mis"][3])
            step = min(step, middle_ref-pos)
        elif error == "ins":
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new += step
        else:
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new -= step

        if error != "ins":
            e_dict[pos] = [error, step]
            pos += step
            if pos >= middle_ref:
                l_new += pos - middle_ref
                middle_ref = pos
        else:
            e_dict[pos - 0.5] = [error, step]

        prev_error = error

        # Randomly select a match length
        for k1 in m_model.keys():
            if k1[0] <= prev_match < k1[1]:
                break
        p = random.random()
        for k2, v2 in m_model[k1].items():
            if k2[0] < p <= k2[1]:
                step = int(np.floor((p - k2[0])/(k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
                break
        # there are no two 0 base matches together
        if prev_match == 0 and step == 0:
            step = 1

        prev_match = step
        if pos + prev_match > middle_ref:
            l_new += pos + prev_match - middle_ref
            middle_ref = pos + prev_match

        pos += prev_match
        if prev_match == 0:
            prev_error += "0"
    return l_new, middle_ref, e_dict

###############################################################################
# FUNCTION mutate_read()    
# INPUTs:
# read: basic read sequence
# read_name
# e_dict
# k: kmer bias  = max length of any homopolymer   
# OUTPUTs:
# read: mutated read sequence  
###############################################################################   
def mutate_read(read, read_name, error_log, e_dict, k, aligned=True):
    search_pattern = "A" * k + "+|" + "T" * k + "+|" + "C" * k + "+|" + "G" * k
    for key in sorted(e_dict.keys(), reverse=True):
        val = e_dict[key]
        key = int(round(key))

        if val[0] == "mis":
            #-----------------------------------------------#
            #-----------------------------------------------#
            ref_base = read[key: key + val[1]]
            while True:
                new_bases = ""
                for i in xrange(val[1]):
                    tmp_bases = list(BASES)
                    tmp_bases.remove(read[key + i])
                    new_base = random.choice(tmp_bases)
                    new_bases += new_base
                check_kmer = read[max(key - k + 1, 0): key] + new_bases + read[key + val[1]: key + val[1] + k - 1]
                if not k or not re.search(search_pattern, check_kmer):
                    break
            new_read = read[:key] + new_bases + read[key + val[1]:]

        elif val[0] == "del":
            new_bases = val[1] * "-"
            ref_base = read[key: key + val[1]]
            new_read = read[: key] + read[key + val[1]:]

        elif val[0] == "ins":
            ref_base = val[1] * "-"
            while True:
                new_bases = ""
                for i in xrange(val[1]):
                    new_base = random.choice(BASES)
                    new_bases += new_base
                check_kmer = read[max(key - k + 1, 0): key] + new_bases + read[key: key + k - 1]
                if not k or not re.search(search_pattern, check_kmer):
                    break
            new_read = read[:key] + new_bases + read[key:]

        read = new_read

        if aligned and val[0] != "match":
            error_log.write(read_name + "\t" + str(key) + "\t" + val[0] + "\t" + str(val[1]) +
                            "\t" + ref_base + "\t" + new_bases + "\n")

    # If choose to have kmer bias, then need to compress homopolymers to 5-mer
    if k:
        read = collapse_homo(read, k)

    return read


def case_convert(seq):
    base_code = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
                 'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T'],
                 'N': ['A', 'T', 'C', 'G'], 'X': ['A', 'T', 'C', 'G']}

    up_string = seq.upper()
    up_list = list(up_string)
    for i in xrange(len(up_list)):
        if up_list[i] in base_code:
            up_list[i] = random.choice(base_code[up_list[i]])
    out_seq = ''.join(up_list)

    return out_seq


def main():
    ref = ""
    model_prefix = "training"
    out = "simulated"
    number = 50
    
    selection_required = False
    # ins, del, mis rate represent the weight tuning in mix model
    kmer_bias = 0
    
    # Parse options and parameters
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hr:c:o:n:i:d:m:", ["select=", "KmerBias=", "seed="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)   
    for opt, arg in opts:
        if opt == "-r":
            ref = arg
        elif opt == "-c":
            model_prefix = arg
        elif opt == "-o":
            out = arg
        elif opt == "-n":
            number = int(arg)
        elif opt == "--select":
            selection_required = True
            selection = arg
        elif opt == "--KmerBias":
            kmer_bias = int(arg)
        elif opt == "--seed":
            random.seed(int(arg))
            np.random.seed(int(arg))
        elif opt == "-h":
            usage()
            sys.exit(0)
        else:
            usage()
            sys.exit(1)
             
    # Generate log file
    #sys.stdout = open(out + ".log", 'w')
    # Record the command typed to log file
    #sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
    #sys.stdout.flush()
        
    print("Reference sequences in file:" + ref)

    if ref == "":
        print("must provide reference sequences!")
        usage()
        sys.exit(1)
             
    nb_seqs = 0
    if selection_required:
        print("Simulation of a specified subset of transcripts \n")
        with open(selection) as in_seqs:
            subset = in_seqs.readlines()
        subset = [x.rstrip('\n') for x in subset] 
        nb_seqs = len(subset)
        sys.stdout.write(strftime("[RNanoSim] %Y-%m-%d %H:%M:%S") + ": Simulation on a subset of " + str(nb_seqs) + " reference sequences. \n")
        #load file and read name by name..?
        #-> save into array of names, check fasta...
        
    # Read in reference genome and generate simulated reads
    print("Reference sequences in file:" + ref)
    print("Number of reads to generate per reference sequence:" + str(number))
    
    sys.stdout.write(strftime("[RNanoSim] %Y-%m-%d %H:%M:%S") + ": Loading model estimates\n")
    trans_read_profile(model_prefix)
    
    # Start simulation
    sys.stdout.write(strftime("[RNanoSim] %Y-%m-%d %H:%M:%S") + ": Start simulation" +"\n")
    sys.stdout.flush()
    out_reads = open(out + "_reads.fasta", 'w')
    
    out_info = open(out + ".read_info.txt", 'w')
    out_info.write("read_name\tread_length\ttrueRef_name\ttrueRef_length\n")

    out_error = open(out + "_error_profile", 'w')
    out_error.write("Seq_name\tSeq_pos\terror_type\terror_length\tref_base\tseq_base\n")
    
    if selection_required:
        out_unfound = open(out + "_unfound.txt", 'w')
    
    seqs_counter = 0
    if selection_required:
        for seq_record in SeqIO.parse(ref, "fasta"):
            s = 0
            while s<len(subset):
                #print(str(s) + " " + str(len(subset)))
                if subset[s]==seq_record.id:
                    #print("Simulating reads for reference sequence:" + seq_record.id)
                    trans_simulation(seq_record, out, kmer_bias, number, out_reads, out_info, out_error)
                    seqs_counter += 1
                    #del subset[s]
                s += 1      
                    
    else:
        for seq_record in SeqIO.parse(ref, "fasta"):
            print("Simulating reads for reference sequence:" + seq_record.id)
            trans_simulation(seq_record, out, kmer_bias, number, out_reads, out_info, out_error)
            seqs_counter += 1
            
    print(str(len(subset)))
    if (len(subset)>0 & selection_required):
        for l in subset:
            out_unfound.write(l+"\n")
        
    in_seqs.close()
    out_unfound.close()
    out_reads.close() 
    out_error.close()
    
    if selection_required:
        sys.stdout.write(strftime("[RNanoSim] %Y-%m-%d %H:%M:%S") + ": "+ str(seqs_counter) + " reference sequences simulated out of " + str(nb_seqs) + " required." +"\n")
    else:    
        sys.stdout.write(strftime("[RNanoSim] %Y-%m-%d %H:%M:%S") + ": "+ str(seqs_counter) + " reference sequences simulated \n")
    sys.stdout.write(strftime("[RNanoSim] %Y-%m-%d %H:%M:%S") + ": Should be finished! \n")
    sys.stdout.close()

if __name__ == "__main__":
    main()
