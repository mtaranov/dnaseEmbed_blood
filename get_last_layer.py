import argparse
import os
import pysam
from pybedtools import BedTool
import numpy as np
from keras import backend as K
from keras.models import Sequential
from keras.layers import Dense
from keras.models import Model
from keras.models import model_from_json
import itertools

def parse_args():
    parser=argparse.ArgumentParser(description='provide json model, weights file, dnase_seq_file, datadir. This script will write last layer embedding  numpy array to datadir')
    parser.add_argument('--json_model',help='json model')
    parser.add_argument('--weights_file',help='weights file')
    parser.add_argument('--dnase_seq_file',help='dnase sequenqe file')
    parser.add_argument('--output_name',help="output name" )
    return parser.parse_args()

def load_model(json_model, weights_file):
    # load json and create model
    json_file = open(json_model, 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    # load weights into new model
    loaded_model.load_weights(weights_file)
    print("Loaded model from disk")
    return loaded_model

def batch_iter(iterable, batch_size):
    '''iterates in batches.
    '''
    it = iter(iterable)
    try:
        while True:
            values = []
            for n in xrange(batch_size):
                values += (it.next(),)
            yield values
    except StopIteration:
        # yield remaining values
        yield values
#        exit 

def generate_from_array(array, batch_size=128):
    """
    Generates the array in batches.
    """
    batch_iterator = batch_iter(array, batch_size)
    for array_batch in batch_iterator:
        yield np.stack(array_batch, axis=0)

def get_last_layer_embedding(data, batch_size,loaded_model):
    generator = generate_from_array(data, batch_size=batch_size)
    intermediate_layer_model = Model(input=loaded_model.layers[0].input,
                                       output=loaded_model.layers[-3].output)
    layer_output = np.vstack([intermediate_layer_model.predict_on_batch(batch) for batch in generator])
    return layer_output


def main():
    args=parse_args()
    X = np.load(args.dnase_seq_file)
    loaded_model = load_model(args.json_model,args. weights_file)
    batch_size = 128
    last_layer = get_last_layer_embedding(X, batch_size,loaded_model)
    
    print "array dim: ", last_layer.shape
    np.save(args.output_name, last_layer)

if __name__=="__main__":
    main()
