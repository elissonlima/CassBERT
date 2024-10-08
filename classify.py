import os
import sys
import subprocess

import pandas as pd

from Bio import SeqIO

def split_file(workdir, inputfile, format, splits):
    fasta_sequences = SeqIO.index(inputfile, format=format)
    keys = list(fasta_sequences.keys())

    chunk_size = int(len(keys) / splits)

    chunks = [[] for _ in range((len(keys) + chunk_size - 1) // chunk_size)]
    for i, item in enumerate(keys):
        chunks[i // chunk_size].append(item)

    result = []

    for i, chunk in enumerate(chunks):
        file_path = os.path.join(workdir, f"tmp_{i}.{format}")
        with open(file_path, "w") as arq:
            for key in chunk:
                SeqIO.write(fasta_sequences[key], arq, format)
            
        result.append(file_path)
    
    return result


if __name__ == '__main__':
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    reference_path = "C:\\Users\\eliss\\Downloads\\AJB1_S86_L007_R1_001.fastq.gz"
    #reference_path = "D:\\code\\CassBERT\\dataset\\eval_2.fasta"
    model_path = "D:\\code\\CassBERT\\trained_models\\protein_final\\model_k5_c3_bp100"
    format = "fastq"
    number_processes = 4
    files = split_file(model_path, reference_path, format, number_processes)
    print(files)
    processes = []
    for idx, file in enumerate(files):
        command = f"d:/code/CassBERT/venv/Scripts/python.exe d:/code/CassBERT/classificator.py {model_path} {file} {idx} {format} 5"

        processes.append(subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True))

    for p in processes:
        p.wait()

    for file in files:
        os.remove(file)

    print("Finished")
        