import logging
import os
import datetime

import numpy as np

import tensorflow as tf

from tensorflow.keras.preprocessing.text import Tokenizer, tokenizer_from_json
from model import Transformer

from Bio import SeqIO
from sklearn.model_selection import train_test_split


class DatasetLoader:

    def __init__(self, k_mer_size=3):
        self.k_mer_size = k_mer_size

    def get_sequences(self, file_obj, size_lenght=34):
        num_reads = 0
        reads = []
        fasta_sequences = SeqIO.parse(file_obj, "fasta")
        for fasta in fasta_sequences:
            sample = []
            # if len(fasta.seq) == size_lenght:
               # and all(c in 'AGCTN' for c in str(fasta.seq).upper()):
            sample.append('<start> ' + ' '.join(str(fasta.seq)[x:x+self.k_mer_size].upper()
                                    for x in range(len(str(fasta.seq))
                                                    - self.k_mer_size + 1)) + ' <end>')
            reads += sample
            print(reads)
            num_reads += len(sample)

        logging.info(f"{num_reads} Reads")

        return reads, num_reads

    def load_data_local(self):
        basedir = os.path.dirname(__file__)
        X_saved_file_name = os.path.join(basedir, 'save', "X_%d.npy" % self.k_mer_size)
        y_saved_file_name = os.path.join(basedir, 'save', "y_%d.npy" % self.k_mer_size)
        tokens_saved_file_name = os.path.join(basedir, 'save', "tokens_%d.txt" % self.k_mer_size)

        X, y = [], []

        if os.path.exists(X_saved_file_name) and \
           os.path.exists( y_saved_file_name):
            X = np.load(X_saved_file_name)
            y = np.load( y_saved_file_name)
        else:       
            arq_cassava = open(os.path.join(basedir, 'dataset','eval_0.fasta'))
            arq_virus = open(os.path.join(basedir, 'dataset','eval_1.fasta'))
            #arq_bacteria = open(os.path.join(basedir, 'dataset','eval_2.fasta'))

            sequences_virus, num_seq_virus = self.get_sequences(arq_virus)
            sequences_cassava, num_seq_cass = self.get_sequences(arq_cassava)
            #sequences_bacteria, num_seq_bact = self.get_sequences(arq_bacteria)

            # labels_virus =    [[0, 1, 0]] * num_seq_virus
            # labels_cassava =  [[1, 0, 0]] * num_seq_cass
            # labels_bacteria = [[0, 0, 1]] * num_seq_bact
            labels_virus =    [[0, 1]] * num_seq_virus
            labels_cassava =  [[1, 0]] * num_seq_cass
            # labels_bacteria = [[0, 0, 1]] * num_seq_bact

            # X = np.array(sequences_virus + sequences_cassava + sequences_bacteria)
            # y = np.array(labels_virus + labels_cassava + labels_bacteria)
            X = np.array(sequences_virus + sequences_cassava)
            y = np.array(labels_virus + labels_cassava)

            arq_virus.close()
            arq_cassava.close()
            np.save(X_saved_file_name, X)
            np.save(y_saved_file_name, y)

        tokenizer = Tokenizer(split=' ')
        if os.path.exists(tokens_saved_file_name):
            with open(tokens_saved_file_name, "r") as arq:
                tokenizer = tokenizer_from_json(arq.read())
        else:
            tokenizer.fit_on_texts(X)

            with open(tokens_saved_file_name, 'w') as arq: 
                arq.write(tokenizer.to_json())
                
        X = np.array(tokenizer.texts_to_sequences(X), dtype=np.object_)
        return X, y, tokenizer


if __name__ == '__main__':
    k_mers = 5
    epochs = 10
    basedir = os.path.dirname(__file__)

    dt = DatasetLoader(k_mer_size=k_mers)
    X, y, tokenizer = dt.load_data_local()

    print("Dataset loaded")
    print("X.shape", X.shape)
    print("y.shape", y.shape)
    batch_size = 64

    maxlen = 102
    vocab_size = len(tokenizer.word_index) + 1
    num_layers = 2
    d_model = 32
    dff = 128
    num_heads = 4
    dropout_rate = 0.1

    model = Transformer(num_layers=num_layers,
                        d_model=d_model,
                        num_heads=num_heads,
                        dff=dff,
                        input_vocab_size=vocab_size,
                        dropout_rate=dropout_rate)
    X_train, X_val, y_train, y_val = train_test_split(X, y,
                                                      test_size=0.33,
                                                      random_state=42)

    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-5),
                  loss=tf.keras.losses.BinaryCrossentropy(),
                  metrics=[tf.keras.metrics.CategoricalAccuracy()])

    current_datetime = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_dir = os.path.join(basedir, 'logs', "k%d-%s" % (int(k_mers), current_datetime))

    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir,
                                                          histogram_freq=1)
    print(X_train.shape)
    print(X_train[0])
    X_train = tf.keras.utils.pad_sequences(X_train, padding="post", maxlen=102, dtype=np.int64)
    X_val = tf.keras.utils.pad_sequences(X_val, padding="post", maxlen=102 , dtype=np.int64)
    print(X_train[0])
    print(X_train.shape)

    print("Starting model training")
    history = model.fit(X_train, y_train,
                        batch_size=batch_size, epochs=epochs,
                        validation_data=(X_val, y_val),
                        verbose=2,
                        callbacks=[tensorboard_callback])
    model.save(os.path.join(basedir, "logs",  "model_k%d_%s" % (int(k_mers), current_datetime)))
    print("Finished")