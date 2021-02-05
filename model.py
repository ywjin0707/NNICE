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
        # subset scdata with selected random column indices (entire batch)
        subset = tf.gather(self.scdata, inpt, batch_dims=1)
        return tf.reduce_sum(subset, axis=1)

class AdversarialSimulator():
    '''
    Ref: github.com/eriklindernoren/Keras-GAN
    '''
    def __init__(self, scdata, n_sim_samples = 1000):
        self.scdata = scdata
        self.n_sim_samples = n_sim_samples
        self.n_celltypes = len(scdata)
        self.n_features = scdata[0].shape[2]
        self.n_cellspersample = 500
        self.n_cellspercelltype = [l.shape[0] for l in scdata]
        optmzr = tf.keras.optimizers.Adam(0.0002, 0.5)
        self.Discriminator = self.build_discriminator()
        self.Discriminator.compile(loss='binary_crossentropy',
            optimizer=optmzr,
            metrics=['accuracy'])
        self.Simulator = self.build_simulator()
        # Simulator takes in indices (from Nprop) as input and generates simbulk
        z = [tf.keras.Input(shape=(500,), dtype=tf.int32, name=f'Input0_{c}') for c in range(self.n_celltypes)]
        simbulk = self.Simulator(z)
        # For the combined model, we only train the Simulator
        self.Discriminator.trainable = False
        # Discriminator takes simbulk as input and determines validity
        valid = self.Discriminator(simbulk)
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
        nprop = nprop * self.n_cellspersample #500
        nprop = tf.cast(nprop, dtype=tf.int32)
        indices = np.zeros((batch_size, 500, self.n_celltypes), dtype=np.int32)
        for c in range(self.n_celltypes):
            for b in range(batch_size):
                # allows for sampling with replacement (increases variability)
                # select {number of cells} of random column indices from scdata with uniform probability
                idx = tf.random.uniform(shape=[nprop.numpy()[b,c]],
                    minval=1, maxval=self.scdata[c].shape[0]-1,
                    dtype=tf.int32)
                indices[b,0:idx.numpy().shape[0],c] = idx
        indices = [indices[:,:,c] for c in range(self.n_celltypes)]
        return nprop, indices
    def build_simulator(self):
        inputs = []
        x = []
        for c in range(self.n_celltypes):
            inpt = tf.keras.Input(shape=(500,), dtype=tf.int32, name=f'Input1_{c}')
            x.append(Subsampling(self.scdata[c])(inpt))
            inputs.append(inpt)
        x = tf.keras.layers.Add()(x)
        # x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Lambda(lambda x: tf.math.log1p(x))(x)
        x = tf.keras.layers.LayerNormalization(center=True, scale=True)(x)
        x = tf.keras.layers.ReLU(max_value=1)(x)
        outputs = tf.keras.layers.Lambda(lambda x: self.MinMaxNorm(x))(x)
        model = tf.keras.Model(inputs, outputs)
        model._name = "Simulator"
        model.summary()
        return model
    def build_discriminator(self):
        inputs = tf.keras.Input(shape=(self.n_features,))
        x = tf.keras.layers.Dense(512)(inputs)
        x = tf.keras.layers.LeakyReLU(alpha=0.2)(x)
        x = tf.keras.layers.Dense(256)(x)
        x = tf.keras.layers.LeakyReLU(alpha=0.2)(x)
        outputs = tf.keras.layers.Dense(1, activation='sigmoid')(x)
        model = tf.keras.Model(inputs, outputs)
        model._name = "Discriminator"
        model.summary()
        return model
    def train(self, X_data, epochs=100, batch_size=32):
        X_data = self.MinMaxNorm(tf.math.log1p(X_data))
        valid = np.ones((batch_size, 1))
        fake = np.zeros((batch_size, 1))
        for epoch in range(epochs):
            '''
            Train Discriminator
            '''
            # Select random subset of X_data equal to batch_size
            idx = np.random.randint(0, X_data.shape[0], batch_size, dtype=np.int64)
            # pdb.set_trace()
            bulk = X_data.numpy()[idx.tolist()]
            # Sample Nprop (cell fractions) using Dirichlet distribution
            nprop, indices = self.simulated_fractions(batch_size)
            # Generate simbulk using Nprop
            simbulk = self.Simulator.predict(indices)
            # Train Discriminator
            d_loss_real = self.Discriminator.train_on_batch(bulk, valid)
            d_loss_fake = self.Discriminator.train_on_batch(simbulk, fake)
            d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)
            '''
            Train Simulator
            '''
            # Train Simulator (wants Discriminator to make mistakes)
            s_loss = self.AdvSim.train_on_batch(indices, valid)
            # Plot the progress
            print ("%d [D loss: %f, acc.: %.2f%%] [G loss: %f]" % (epoch, d_loss[0], 100*d_loss[1], s_loss))
