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

def vae_model(input_size=11084, intermediate_size=2000, latent_size=150):
    """
    Setup model
    """
    # VAE model = encoder + decoder
    # build encoder model
    inputs = Input(shape=(input_size,), name='encoder_input')
    x = Dense(intermediate_size, activation='relu')(inputs)
    x = Dense(intermediate_size, activation='relu')(x)
    z_mean = Dense(latent_size, name='z_mean')(x)
    z_log_var = Dense(latent_size, name='z_log_var')(x)
    # use reparameterization trick to push the sampling out as input
    z = Lambda(sampling, name='z')([z_mean, z_log_var])
    # instantiate encoder model
    encoder = Model(inputs, [z_mean, z_log_var, z], name='encoder')
    # build decoder model
    latent_inputs = Input(shape=(latent_size,), name='z_sampling')
    x = Dense(intermediate_size, activation='relu')(latent_inputs)
    x = Dense(intermediate_size, activation='relu')(x)
    outputs = Dense(input_size, activation='sigmoid')(x)
    decoder = Model(latent_inputs, outputs, name='decoder')
    # instantiate VAE model
    outputs = decoder(encoder(inputs)[2])
    # Loss
    reconstruction_loss = mse(inputs, outputs) * input_size
    kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    vae_loss = K.mean(reconstruction_loss + kl_loss)
    vae = Model(inputs, outputs, name='vae')
    vae.add_loss(vae_loss)
    return vae, encoder, decoder

def format_preds(pred_mean, pred_log_var, strains):
    """
    Format model predictions into a pandas DataFrame
    """
    out = np.empty((pred_mean.shape[0], pred_mean.shape[1] + pred_log_var.shape[1]),
                   dtype=pred_mean.dtype)
    out[:, 0::2] = pred_mean
    out[:, 1::2] = pred_log_var
    colnames = [f'{j}_{i}' for i in range(1, pred_mean.shape[1] + 1)
                for j in ('mean', 'log_var')]
    df = pd.DataFrame(data=out, columns=colnames)
    df['strain'] = strains
    df = df[['strain'] + colnames]
    return df

def main(args):
    """Main"""
    train = pd.read_csv(args.train, sep='\t', index_col=0)
    test = pd.read_csv(args.test, sep='\t', index_col=0)
    train_mat = train.to_numpy().astype(np.float32)
    test_mat = test.to_numpy().astype(np.float32)

    callbacks = []
    callbacks.append(keras.callbacks.TensorBoard(log_dir=args.outdir, histogram_freq=50,
                                                 write_graph=True, profile_batch=0))

    vae, encoder, decoder = vae_model()
    vae.compile(optimizer='adam')
    vae.fit(train_mat, epochs=args.epochs, batch_size=100,
            validation_data=(test_mat, None), callbacks=callbacks)

    train_pred = encoder(train_mat)[:2]
    test_pred = encoder(test_mat)[:2]

    train_pred = format_preds(train_pred[0].numpy(), train_pred[1].numpy(), list(train.index))
    test_pred = format_preds(test_pred[0].numpy(), test_pred[1].numpy(), list(test.index))
    preds = pd.concat([train_pred, test_pred])

    vae.save(f'{args.outdir}/vae.tf', save_format='tf')
    encoder.save(f'{args.outdir}/encoder.tf', save_format='tf')
    preds.to_csv(f'{args.outdir}/profiles.tsv', index=False, sep='\t')

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--train', '-r', default='data/full_omics_train.tsv',
                        help="Training omics data file")

    parser.add_argument('--test', '-e', default='data/full_omics_test.tsv',
                        help="Testing omics data file")

    parser.add_argument('--epochs', '-p', default=100, type=int, help="Epochs to train for")

    parser.add_argument('--outdir', '-o', default='data/vae', help="Output directory")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())