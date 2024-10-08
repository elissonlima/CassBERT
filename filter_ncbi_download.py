import os
import sys
import random

from Bio import SeqIO


if __name__ == '__main__':
    path = 'D:\\docker_projects\\bioinformatics\\code\\bacterias'
    output_file = open("D:\\code\\CassBERT\\references\\bacterias_protein.fasta", "a")
    for subdir, dirs, files in os.walk(path):
        for file in files:
            file_name = os.path.join(subdir, file)
            if file_name.endswith("faa"):
                with open(file_name) as arq:
                    fasta_sequences = SeqIO.index(file_name, "fasta")
                    ids_ = random.choices(list(fasta_sequences.keys()), k=50)
                    for id_ in ids_:
                        SeqIO.write(fasta_sequences[id_], output_file, "fasta")
            else:
                print(file_name)
                    