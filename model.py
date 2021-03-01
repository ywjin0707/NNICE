import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow.keras import backend as K
import numpy as np
import pdb

class Subsampling(tf.keras.layers.Layer):

    def __init__(self, scdata):
        super(Subsampling, self).__init__()
        # initialize one layer for each cell type
        self.scdata=scdata

    def call(self, inpt):
        # select {number of cells} of random column indices from scdata with uniform probability
        # allows for sampling with replacement (increases variability)
        idx = tf.random.uniform(
            shape=tf.squeeze(inpt),
            minval=0, maxval=self.scdata.shape[0],
            dtype=tf.int32)
        # subset scdata with selected random column indices
        subset = tf.gather(self.scdata, idx, axis=0)
        return tf.reduce_sum(subset, axis=0)

class AdversarialSimulator():

    '''
    Ref: github.com/eriklindernoren/Keras-GAN
    '''

    def __init__(self, scdata, n_sim_samples = 1000):
        self.scdata = scdata
        self.n_sim_samples = n_sim_samples
        self.n_celltypes = len(scdata)
        self.n_features = scdata[0].shape[1]
        self.n_cellspersample = 500
        self.n_cellspercelltype = [l.shape[0] for l in self.scdata]
        optmzr = tf.keras.optimizers.Adam(0.0002, 0.5)
        self.Discriminator = self.build_discriminator()
        self.Discriminator.compile(loss='binary_crossentropy',
            optimizer=optmzr,
            metrics=['accuracy'])
        self.Simulator = self.build_simulator()
        # Simulator takes in Nprop as input and generates simbulk
        z = tf.keras.layers.Input(shape=(self.n_celltypes), dtype=tf.int32)
        img = self.Simulator(z)
        # For the combined model, we only train the Simulator
        self.Discriminator.trainable = False
        # Discriminator takes simbulk as input and determines validity
        valid = self.Discriminator(self.simbulk)
        # stacked Simulator and Discriminator
        self.AdvSim = tf.keras.Model(z, valid)
        self.AdvSim.compile(loss='binary_crossentropy', optimizer=optmzr)

    def MinMaxNorm(self, x):
        x_scaled = tf.math.divide_no_nan(
            (x - tf.math.reduce_min(x)),
            (tf.math.reduce_max(x) - tf.math.reduce_min(x)))
        return x_scaled

    def simulated_fractions(self, batch_size):
        alpha = [1]*self.n_celltypes
        dist = tfp.distributions.Dirichlet(alpha)
        nprop = dist.sample([batch_size])
        nprop = nprop * self.n_cellspersample
        nprop = tf.cast(nprop, dtype=tf.int32)
        return nprop

    def build_simulator(self):
        inputs = []
        x = []
        for c in range(self.n_celltypes):
            inpt = tf.keras.layers.Input(shape=(1), dtype = tf.int32, name=f'input{c}')
            x.append(Subsampling(self.scdata[c])(inpt))
            inputs.append(inpt)
        x = tf.keras.layers.Add()(x)
        x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Lambda(lambda x: tf.math.log1p)(x)
        x = tf.keras.layers.LayerNormalization(center=True, scale=True)(x)
        outputs = tf.keras.layers.Lambda(lambda x: self.MinMaxNorm)(x)
        model = tf.keras.Model(inputs, outputs)
        model.summary()
        return model

    def build_discriminator(self):
        inputs = tf.keras.Input(shape=(self.n_features,))
        x = tf.keras.layers.Dense(256)(inputs)
        x = tf.keras.layers.LeakyReLU(alpha=0.2)(x)
        x = tf.keras.layers.Dense(128)(x)
        x = tf.keras.layers.LeakyReLU(alpha=0.2)(x)
        outputs = tf.keras.layers.Dense(1, activation='sigmoid')(x)
        model = tf.keras.Model(inputs, outputs)
        model._name = "Discriminator"
        model.summary()
        return model

    def train(self, X_data, steps=1000):
        X_data = self.MinMaxNorm(tf.math.log1p(X_data))
        valid = np.ones((batch_size, 1))
        fake = np.zeros((batch_size, 1))
        for step in range(steps):
            '''
            Train Discriminator
            '''
            # Select random subset of X_data equal to batch_size
            idx = np.random.randint(0, X_data.shape[0], 1)
            bulk = X_data.numpy()[idx]
            # Sample Nprop (cell fractions) using Dirichlet distribution
            nprop = self.simulated_fractions(batch_size)
            # Generate simbulk using Nprop
            simbulk = self.Simulator.predict(nprop)
            # Train Discriminator
            d_loss_real = self.Discriminator.train_on_batch(bulk, valid)
            d_loss_fake = self.Discriminator.train_on_batch(simbulk, fake)
            d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)
            '''
            Train Simulator
            '''
            # Train Simulator (wants Discriminator to make mistakes)
            s_loss = self.AdvSim.train_on_batch(nprop, valid)
            # Plot the progress
            print ("Step: %d [D loss: %f, acc.: %.2f%%] [G loss: %f]" % (step, d_loss[0], 100*d_loss[1], s_loss))
