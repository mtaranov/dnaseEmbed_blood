import argparse
import os
import pysam
from pybedtools import BedTool
import numpy as np

def parse_args():
    parser=argparse.ArgumentParser(description='provide dnase_narrowPeak_file, source fasta file, window size around peak summit, datadir. This script will write seq_at_dnasePeaks numpy array to datadir')
    parser.add_argument('--fasta_file',help='reference fasta file')
    parser.add_argument('--dnase_narrowPeak_file',help='dnase_narrowPeak file')
    parser.add_argument('--window',type=int,help='window size aeound dnase_narrowPeak summit')
    parser.add_argument('--output_name',help="output name" )
    return parser.parse_args()

def seq_to_one_hot(sequence): 
    #assumes 1d and tensorflow dim ordering
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(to_return, sequence, one_hot_axis=1)
    return to_return

# Letter as 1, other letters as 0
def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #zeros_array should be an array of dim 4xlen(sequence), filled with zeros.
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1

def get_seq_to_one_hot_at_window_around_dnasePeaks(dnase_narrowPeak_file, fasta_file, window):
    data = BedTool(dnase_narrowPeak_file)
    count = 0
    for i in data:
        chrom = str(i[0])
        pos_start, pos_end = int(i[1])+int(i[9])-window/2, int(i[1])+int(i[9])+window/2
        seq = pysam.FastaFile(fasta_file).fetch(chrom, pos_start, pos_end)
        seq_array_for_line = seq_to_one_hot(seq)
        if count == 0:
            seq_array = seq_array_for_line[None,:,:]
        else: 
            seq_array = np.concatenate((seq_array,seq_array_for_line[None,:,:]), axis=0)
        count += 1
    return seq_array

# Pyton's 
def loadSequencesfromFasta(inputFasta, indexForAnalysis=None):
    t0 = time.time()
    if indexForAnalysis is not None: 
        sequences = [el.strip('\n') for (i, el) in enumerate(os.popen('zcat {0} | grep -v ">"'.format(
            inputFasta)).readlines()) if i in indexForAnalysis]
    else:
        sequences = [el.strip('\n') for el in os.popen('zcat {0} | grep -v ">"'.format(
            inputFasta)).readlines()]
    inputSeq = util.setOfSeqsTo2Dimages(sequences)
    inputSeq = inputSeq.astype('float32')
    print 'Loaded Input Sequences'
    t1 = time.time()
    print 'Time: %s'%(t1-t0)
    return inputSeq

def main():
    args=parse_args()
    seq_array = get_seq_to_one_hot_at_window_around_dnasePeaks(args.dnase_narrowPeak_file, args.fasta_file, args.window)
    print "array dim: ", seq_array.shape
    np.save(args.output_name, seq_array)

if __name__=="__main__":
    main()
