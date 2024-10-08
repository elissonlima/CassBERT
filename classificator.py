import os
import tensorflow as tf
import numpy as np
from tensorflow.keras.preprocessing.text import Tokenizer, tokenizer_from_json
from Bio import SeqIO
import time
from datetime import  timedelta
import sys
import logging
import pandas as pd

from biogl import translate, rev_comp
from collections import defaultdict


class NucleotideClassificator:

    def __init__(self, model_path, k_mers_size):
        self.k_mer_size = k_mers_size
        self.model = tf.keras.models.load_model(os.path.join(model_path, 'model'), compile=False)
        #self.model = tf.keras.layers.TFSMLayer(os.path.join(model_path, 'model'), call_endpoint="serving_default")
        with open(os.path.join(model_path, "tokens.txt"), "r") as arq:
            self.tokenizer = tokenizer_from_json(arq.read())
    

    def _input_to_token(self, input_value):
        #kmers_input = '<start> '.join(input_value[x:x+self.k_mer_size].upper() for x in range(len(input_value) - self.k_mer_size + 1)) + ' <end>'

        kmers_input = "<start> "
        kmers_input += ' '.join([input_value[x:x+self.k_mer_size].upper() for x in range(len(input_value) - self.k_mer_size + 1)][:100])
        kmers_input += " <end>"
        
        # kmers_input =[input_value[x:x+self.k_mer_size].upper() for x in range(len(input_value) - self.k_mer_size + 1)]
        # if self.k_mer_size == 3:
        #     return tf.constant(self.tokenizer.texts_to_sequences([kmers_input]), dtype=tf.float32)
        # else: 
        return tf.constant(self.tokenizer.texts_to_sequences([kmers_input]), dtype=tf.int64)
    def __call__(self, input_value):
        try:
            tokenized_input = self._input_to_token(input_value)
            #padded = tf.keras.utils.pad_sequences(tokenized_input, maxlen=30, value=0, dtype='int64')
            padded = tf.keras.utils.pad_sequences(tokenized_input, padding="post", maxlen=102, dtype=np.int64)
            
            res = self.model(padded).numpy()
            # if np.argmax(res[0]) == 1:
            #     return 'virus'
            # elif np.argmax(res[0]) == 0:
            #     return 'cassava'
            return res[0]
        except Exception as e:
            print(e)
            logging.info(e)
        

def run_model(model_path, k_mers_size, reference_path, process_idx, format, sequence_max_size=34):    
    logging.info("Starting for file %s" % reference_path)
    filename=os.path.basename(reference_path)

    model = NucleotideClassificator(model_path, k_mers_size)
    fasta_sequences = SeqIO.index(reference_path, format=format)
    keys = list(fasta_sequences.keys())

    total_amount = len(keys)
    count_eval = 0
    start_time = time.time()

    # arq_virus = 
    # arq_cassava = open(os.path.join(model_path, f"res_cassava.fasta"), "w")

    dict_virus = {'seq_id':[], 'prob':[]}
    num_virus = 0
    dict_bacteria = {'seq_id':[], 'prob':[]}
    num_bacteria = 0
    dict_cassava = {'seq_id':[], 'prob':[]}
    num_cassava = 0
    
    for idx in keys:

        sequence = str(fasta_sequences[idx].seq)
        translated_0 = translate(sequence, phase=0)
        translated_1 = translate(sequence, phase=1)
        translated_2 = translate(sequence, phase=2)
        rev_translated_0 =  translate(rev_comp(sequence), phase=0)
        rev_translated_1 =  translate(rev_comp(sequence), phase=1)
        rev_translated_2 =  translate(rev_comp(sequence), phase=2)
        has_virus = False
        virus_probs = [-1]
        has_bacteria = False
        bacteria_probs = [-1]
        cassava_probs = [-1]

        for protein_seq in [translated_0, translated_1, translated_2,
                            rev_translated_0, rev_translated_1,
                            rev_translated_2]:
            
            if '*' in protein_seq:
                continue
            
            res = model(protein_seq)
            # print(res)

            if np.argmax(np.array(res)) == 1:
                has_virus = True
                virus_probs.append(res[1])
            elif np.argmax(np.array(res)) == 2:
                has_bacteria = True
                bacteria_probs.append(res[2])
            else:
                cassava_probs.append(res[0])

        # if has_virus:
        #     dict_virus['seq_id'].append(fasta_sequences[idx].id)
        #     dict_virus['prob'].append(max(virus_probs))
        #     num_virus += 1
        # else:
        #     dict_cassava['seq_id'].append(fasta_sequences[idx].id)
        #     dict_cassava['prob'].append(max(cassava_probs))
        #     num_cassava += 1

        if has_virus and not has_bacteria:
            #dict_virus[fasta_sequences[idx].id] = max(virus_probs)
            dict_virus['seq_id'].append(fasta_sequences[idx].id)
            dict_virus['prob'].append(max(virus_probs))
            num_virus += 1
        elif has_bacteria and not has_virus:
            dict_bacteria['seq_id'].append(fasta_sequences[idx].id)
            dict_bacteria['prob'].append(max(bacteria_probs))
            num_bacteria +=1
        elif has_virus and has_bacteria:
            if max(virus_probs) > max(bacteria_probs):
                dict_virus['seq_id'].append(fasta_sequences[idx].id)
                dict_virus['prob'].append(max(virus_probs))
                num_virus += 1 
            else:
                dict_bacteria['seq_id'].append(fasta_sequences[idx].id)
                dict_bacteria['prob'].append(max(bacteria_probs))
                num_bacteria += 1 
        else:
            dict_cassava['seq_id'].append(fasta_sequences[idx].id)
            dict_cassava['prob'].append(max(cassava_probs))
            num_cassava += 1

        count_eval += 1

        if (count_eval % 100000) == 0:
            elapsed = time.time() - start_time
            logging.info("{} - {}/{} Sequences evaluated - Virus: {} | Bacteria: {} | Cassava: {} - Elapsed: {}".format(filename, count_eval, total_amount, num_virus, num_bacteria, num_cassava, timedelta(seconds=elapsed)))
            start_time = time.time()

    elapsed = time.time() - start_time
    logging.info("{} - {}/{} Sequences evaluated - Virus: {} | Bacteria: {} | Cassava: {} - Elapsed: {}".format(filename, count_eval, total_amount, num_virus, num_bacteria, num_cassava, timedelta(seconds=elapsed)))

    df_virus  = pd.DataFrame.from_dict(dict_virus)
    df_virus.to_csv(os.path.join(model_path, f"virus_{process_idx}.csv"), index=False, header=False)

    df_cassava  = pd.DataFrame.from_dict(dict_cassava)
    df_cassava.to_csv(os.path.join(model_path, f"cassava_{process_idx}.csv"), index=False, header=False)

    df_bacteria  = pd.DataFrame.from_dict(dict_bacteria)
    df_bacteria.to_csv(os.path.join(model_path, f"bacteria_{process_idx}.csv"), index=False, header=False)

if __name__ == '__main__':
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

    model_path = sys.argv[1]
    input_file_path = sys.argv[2]
    process_idx = sys.argv[3]
    format = sys.argv[4]
    k_mers_size = int(sys.argv[5])


    logging.basicConfig(filename=f"process_{process_idx}_log.out",
                filemode='w',
                format='%(asctime)s - %(message)s',
                datefmt='%H:%M:%S',
                level=logging.DEBUG)
    run_model(model_path, k_mers_size, input_file_path, process_idx, format)
    
