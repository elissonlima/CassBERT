# Create traninig dataset from "references" folder

import os
import pathlib
from operator import itemgetter
from typing import Sequence
from Bio import SeqIO
from Bio.Seq import Seq
import random
import numpy as np
import logging
import pathlib
from biogl import rev_comp

#from model.util import timer

# Creating constants

OUT_SEQUENCES_LENTH = 100
DEPTH = "max"

REFERENCE_PATH = 'D:\\code\\cassava-virus-finder\\data\\references_virus'

# Function that generate the train and test datasets containg genome random
# reads
#@timer
def generate_datasets_virus(read_size=100):
    
    read_count = 0

    for file_name in os.listdir(REFERENCE_PATH):

        destination_path = "D:\\code\\cassava-virus-finder\\data\\"\
            "train_%d\\virus\\%s" % (read_size, file_name)
        pathlib.Path(os.path.dirname(destination_path)).mkdir(parents=True, exist_ok=True)
        
        fpath = os.path.join(REFERENCE_PATH, file_name)

        cnt = generate_reads_using_sliding_window(fpath, read_size, destination_path,
                                            incluse_reverse_complement=True)

        logging.info("file: %s [ OK ] - %d Reads generated" % (file_name, cnt))
        read_count += cnt

    logging.info("%d Reads generated for Virus class" % read_count)

#@timer
def generate_datasets_cassava(read_size=100):
    
    # Generate cassava reads
    cassava_file = 'Mesculenta_671_v8.0.fna'
    fpath = 'D:\\code\\cassava-virus-finder\\data\\references_cassava\\'\
    'Mesculenta_671_v8.0.fna'

    destination_path = f"D:\\code\\cassava-virus-finder\\data\\"\
                "train_%d\\cassava\\" % (read_size)
    pathlib.Path(destination_path).mkdir(parents=True, exist_ok=True)

    cnt = generate_reads_using_sliding_window_split_files(fpath, read_size, destination_path,
                                            incluse_reverse_complement=False)

    logging.info("file: %s [ OK ] - %d Reads generated" % (cassava_file, cnt))
    logging.info("%d Reads generated for Cassava class" % cnt)


def generate_reads_using_sliding_window(in_file, out_len, out_file,
                                        incluse_reverse_complement=True):
    

    fasta_sequences = SeqIO.index(in_file, "fasta")
    read_count = 0

    with open(out_file, "w") as arq:
        for idx in fasta_sequences.keys():

            sequence_obj = fasta_sequences[idx]
            sequence = str(sequence_obj.seq)
            seq_name = sequence_obj.id

            for i in range(len(sequence) - (out_len - 1)):
                read = sequence[i:i+out_len]
                read_id = "%s:%d:%d" % (seq_name,i,-1)

                # Check if read isn't empty
                if read:
                    arq.write(">%s\n%s\n" % (read_id, read))
                    read_count += 1

                    if incluse_reverse_complement:
                        seq = Seq(read)
                        reversed_read = str(seq.reverse_complement())
                        reversed_read_id = "%s_reversed:%d:%d" % (seq_name,i,-1)
                        arq.write(">%s\n%s\n" % (reversed_read_id, reversed_read))
                        read_count += 1
    return read_count


def generate_reads_using_sliding_window_split_files(in_file, out_len, out_path,
                                        incluse_reverse_complement=True):
    

    fasta_sequences = SeqIO.index(in_file, "fasta")
    read_count = 0
    file_idx = 0
    arq = open(os.path.join(out_path, "file_%d.fasta" % file_idx), "a")

    for idx in fasta_sequences.keys():

        sequence_obj = fasta_sequences[idx]
        sequence = str(sequence_obj.seq)
        seq_name = sequence_obj.id        

        for i in range(0, len(sequence) - (out_len - 1), 100):
            read = sequence[i:i+out_len]
            read_id = "%s:%d:%d" % (seq_name,i,-1)
            if read_count % 10e7 == 0:
                file_idx+=1
                arq.close()
                arq = open(os.path.join(out_path, "file_%d.fasta" % file_idx), "a")
            # Check if read isn't empty
            if read:
                arq.write(">%s\n%s\n" % (read_id, read))
                read_count += 1

                # if incluse_reverse_complement:
                #     seq = Seq(read)
                #     reversed_read = str(seq.reverse_complement())
                #     reversed_read_id = "%s_reversed:%d:%d" % (seq_name,i,-1)
                #     arq.write(">%s\n%s\n" % (reversed_read_id, reversed_read))
                #     read_count += 1
    arq.close()
    return read_count


# # Read file and generate random reads
# def generate_reads_from_file(file_path, read_size, file_type="fasta",
#                              generate_type="window"):
#     '''
#     Read file and generate random reads

#         Parameters:
#             file_path (str): Path to file to get reads
#             file_type (str): File type
#         Returns:
#     '''   
#     fasta_sequences = SeqIO.parse(open(file_path),file_type)
#     result_reads = []
#     for read in fasta_sequences:
#         if generate_type == "window":
#             result_reads += generate_reads_using_sliding_window(str(read.seq), 
#                                                                 read.id,
#                                                                 read_size)
#         elif generate_type == "random":            
#             result_reads += get_random_reads_from_seq(str(read.seq), read.id, 
#                                 read_size, 1)

#     return result_reads


def get_random_reads_from_seq(input_file, output_file, max_read_len, depth, strand='both'):
    '''
    Generate random reads from sequence string

        Parameters:
            sequence (str): Sequence as string
            seq_name (str): Sequence id
            max_read_len (str): Max size of output reads
            depth:
                int value -> How many time Increase the number of out reads
                "max" -> Get maximum distinct reads from sequence
            strand (str): 
                both -> generate forward or reverse read (randomly)
                antisense -> only generate reverse read
        Returns:
    '''  

    fasta_sequences = SeqIO.index(input_file, "fasta")
    arq_out =  open(output_file, "w")

    for idx in fasta_sequences.keys():
        
        sequence = str(fasta_sequences[idx].seq)

        maxstart = len(sequence) - max_read_len + 1
        maxbp = len(sequence) # max number of bp to be generated
        if depth == "max":
            maxbp = float('inf')
        else:
            maxbp *= depth
        #$thisseq =~ tr/ACGT/TGCA/;
        bpcount = 0
        repetead_count = 0

        result_set = set()

        if len(sequence) < max_read_len:
            continue
        elif len(sequence) == max_read_len:
            arq_out.write(">%s\n%s\n" % (fasta_sequences[idx].id, 
                                         sequence))
        
        while bpcount < maxbp and repetead_count < len(sequence):

            startidx = random.randint(0, maxstart)
            new_seq = sequence[startidx:startidx+max_read_len]

            num_strand = 1
            if strand == "both":
                if startidx % 2 == 0:
                    num_strand = -1
                    new_seq = rev_comp(new_seq)
                    # new_seq = new_seq[::-1]                    
                    # translate_dict = new_seq.maketrans("ACGT","TGCA")
                    # new_seq =   new_seq.translate(translate_dict)
            elif strand == "antisense":
                new_seq = rev_comp(new_seq)
            
                # new_seq = new_seq[::-1]
                # translate_dict = new_seq.maketrans("ACGT","TGCA")
                # new_seq = new_seq.translate(translate_dict)

            if new_seq in result_set:
                repetead_count += 1
            else:
                repetead_count = 0
                result_id = "%s:%d:%d" % (fasta_sequences[idx].id,startidx,num_strand)
                arq_out.write(">%s\n%s\n" % (result_id, new_seq))

            bpcount += max_read_len

def get_random_reads_from_seq_dif_sizes(input_file, output_file, sizes, depth, strand='both'):
    '''
    Generate random reads from sequence string

        Parameters:
            sequence (str): Sequence as string
            seq_name (str): Sequence id
            max_read_len (str): Max size of output reads
            depth:
                int value -> How many time Increase the number of out reads
                "max" -> Get maximum distinct reads from sequence
            strand (str): 
                both -> generate forward or reverse read (randomly)
                antisense -> only generate reverse read
        Returns:
    '''  

    fasta_sequences = SeqIO.index(input_file, "fasta")
    arq_out =  open(output_file, "w")

    for idx in fasta_sequences.keys():

        result_set = set()
        repetead_count = 0

        for max_read_len in sizes:

            bpcount = 0
        
            sequence = str(fasta_sequences[idx].seq)

            maxstart = len(sequence) - max_read_len + 1
            maxbp = len(sequence) # max number of bp to be generated
            if depth == "max":
                maxbp = float('inf')
            else:
                maxbp *= depth
            #$thisseq =~ tr/ACGT/TGCA/;
        
            if len(sequence) < max_read_len:
                continue
            elif len(sequence) == max_read_len:
                arq_out.write(">%s\n%s\n" % (fasta_sequences[idx].id, 
                                            sequence))
                continue
            
            while bpcount < maxbp and repetead_count < len(sequence):

                startidx = random.randint(0, maxstart)
                new_seq = sequence[startidx:startidx+max_read_len]

                num_strand = 1
                # if strand == "both":
                #     if startidx % 2 == 0:
                #         num_strand = -1
                #         new_seq = rev_comp(new_seq)
                #         # new_seq = new_seq[::-1]                    
                #         # translate_dict = new_seq.maketrans("ACGT","TGCA")
                #         # new_seq =   new_seq.translate(translate_dict)
                # elif strand == "antisense":
                #     new_seq = rev_comp(new_seq)
                    # new_seq = new_seq[::-1]
                    # translate_dict = new_seq.maketrans("ACGT","TGCA")
                    # new_seq = new_seq.translate(translate_dict)

                #if new_seq in result_set:
                if any(new_seq in s for s in result_set):
                    repetead_count += 1
                else:
                    repetead_count = 0
                    result_set.add(new_seq)
                    result_id = "%s:%d:%d" % (f"{fasta_sequences[idx].id}_{max_read_len}",startidx,num_strand)
                    arq_out.write(">%s\n%s\n" % (result_id, new_seq))

                bpcount += max_read_len

def remove_duplicates_from_file(inputfile, output_file):

    sequences = SeqIO.parse(input_file, "fasta")
    arq_out = open(output_file, "w")
    id_set = set()
    for seq in sequences:
        if seq.id not in id_set:
            SeqIO.write(seq, arq_out, "fasta")
            id_set.add(seq.id)

        
        
if __name__ == '__main__':
    # # cassava depth 60
    # # virus depth 30
    # # bacteria depth 2
    input_file = "D:\\code\\cassava-virus-finder\\data\\references_cassava\\Mesculenta_671_v8.0.fna"
    output_file = "D:\\code\\CassBERT\\dataset\\cassava_nucleotide.fasta"
    get_random_reads_from_seq(input_file, output_file, 100, 1, strand=None)
    #get_random_reads_from_seq_dif_sizes(input_file, output_file, [34, 50, 100], 1, strand=None)
    # # output_file_no_duplicates = "D:\\code\\CassBERT\\references\\cassava_protein.fasta"
    # output_file = "D:\\code\\CassBERT\\references\\bacteria_nodup_protein.fasta"
    # get_random_reads_from_seq_dif_sizes(input_file, 
    #                           output_file, 
    #                           [34, 50, 100],
    #                           1, strand='both')
    #remove_duplicates_from_file(input_file, output_file)
    


