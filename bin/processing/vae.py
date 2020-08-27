#!/usr/bin/env python3
"""
Run a VAE on the omics data to generate latent representations

VAE adapted from Keras examples (https://github.com/keras-team/keras/blob/master/examples/variational_autoencoder.py)
"""
import argparse
import numpy as np
import pandas as pd

import tensorflow.keras as keras
from tensorflow.keras.layers import Lambda, Input, Dense
from tensorflow.keras.models import Model
from tensorflow.keras.losses import mse, binary_crossentropy
from tensorflow.keras import backend as K

# reparameterization trick
# instead of sampling from Q(z|X), sample epsilon = N(0,I)
# z = z_mean + sqrt(var) * epsilon
def sampling(args):
    """Reparameterization trick by sampling from an isotropic unit Gaussian.
    # Arguments
        args (tensor): mean and log of variance of Q(z|X)
    # Returns
        z (tensor): sampled latent vector
    """
    z_mean, z_log_var = args
    batch = K.shape(z_mean)[0]
    dim = K.int_shape(z_mean)[1]
    # by default, random_normal has mean = 0 and std = 1.0
    epsilon = K.random_normal(shape=(batch, dim))
    return z_mean + K.exp(0.5 * z_log_var) * epsilon

def vae_model(input_size=1788, intermediate_size=500, latent_size=100):
    """
    Setup model
    """
    # VAE model = encoder + decoder
    # build encoder model
    inputs = Input(shape=(input_size,), name='encoder_input')
    x = Dense(intermediate_size, activation='relu')(inputs)
    z_mean = Dense(latent_size, name='z_mean')(x)
    z_log_var = Dense(latent_size, name='z_log_var')(x)
    # use reparameterization trick to push the sampling out as input
    z = Lambda(sampling, name='z')([z_mean, z_log_var])
    # instantiate encoder model
    encoder = Model(inputs, [z_mean, z_log_var, z], name='encoder')
    # build decoder model
    latent_inputs = Input(shape=(latent_size,), name='z_sampling')
    x = Dense(intermediate_size, activation='relu')(latent_inputs)
    outputs = Dense(input_size, activation='sigmoid')(x)
    decoder = Model(latent_inputs, outputs, name='decoder')
    # instantiate VAE model
    outputs = decoder(encoder(inputs)[2])
    return Model(inputs, outputs, name='vae')

def main(args):
    """Main"""
    data = pd.read_csv(args.omics, sep='\t')
    model = vae_model()


def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--omics', '-o', default='data/full_omics.tsv', help="Omics data file")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())