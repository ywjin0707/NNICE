import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow.keras import backend as K

class Subsampling(tf.keras.layers.Layer):
    def __init__(self, scdata):
        super(Subsampling, self).__init__()
        # initialize one layer for each cell type
        self.scdata=scdata
    def call(self, inputs):
        # select {number of cells} of random column indices from scdata with uniform probability
        # allows for sampling with replacement (increases variability)
        idx = tf.random.uniform(
            shape=[inputs],
            minval=0, maxval=self.scdata.shape[1]-1,
            dtype=tf.int32)
        # subset scdata with selected random column indices
        subset = tf.gather(self.scdata, idx, axis=1)
        return tf.reduce_sum(subset, axis=1)

class AdversarialSimulator():
    '''
    Ref: github.com/eriklindernoren/Keras-GAN
    '''
    def __init__(self, scdata, n_sim_samples = 1000):
        self.scdata = scdata
        self.n_sim_samples = n_sim_samples
        self.n_celltypes = len(scdata)
        self.n_features = scdata[0].shape[1]
        optmzr = tf.keras.optimizers.Adam(0.0002, 0.5)
        self.Discriminator = self.build_discriminator()
        self.Discriminator.compile(loss='binary_crossentropy',
            optimizer=optmzr,
            metrics=['accuracy'])
        self.Simulator = self.build_simulator()
        # Simulator takes in Nprop as input and generates simbulk
        z = Input(shape=(self.n_celltypes))
        img = self.Simulator(z)
        # For the combined model, we only train the Simulator
        self.Discriminator.trainable = False
        # Discriminator takes simbulk as input and determines validity
        valid = self.Discriminator(self.simbulk)
        # stacked Simulator and Discriminator
        self.AdvSim = tf.keras.Model(z, valid)
        self.AdvSim.compile(loss='binary_crossentropy', optimizer=optmzr)

    def MinMaxNorm(x):
        x_scaled = tf.math.divide_no_nan(
            (x - tf.math.reduce_min(x)),
            (tf.math.reduce_max(x) - tf.math.reduce_min(x)))
        return x_scaled

    def simulated_fractions(self, batch_size):
        alpha = [1]*self.n_celltypes
        dist = tfp.distributions.Dirichlet(alpha)
        nprop = dist.sample([batch_size])
        return nprop

    def build_simulator(self):
        inputs = tf.keras.layers.Input(shape=self.n_celltypes)
        x = []
        for c in range(n_celltype):
            x.append(Subsampling(self.scdata[c])(inputs))
        x = tf.keras.layers.Add()(x)
        x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Lambda(lambda x: tf.math.log1p)(x)
        x = tf.keras.layers.LayerNormalization(center=True, scale=True)(x)
        outputs = tf.keras.layers.Lambda(lambda x: self.MinMaxNorm)(x)
        model = tf.keras.Model(inputs, outputs)
        model.summary()
        return model

    def build_discriminator(self):
        inputs = tf.keras.layers.Input(shape=self.n_features)
        x = tf.keras.layers.Dense(512)(inputs)
        x = tf.keras.layers.LeakyReLU(alpha=0.2)(x)
        x = tf.keras.layers.Dense(256)(x)
        x = tf.keras.layers.LeakyReLU(alpha=0.2)(x)
        x = tf.keras.layers.Dense(1, activation='sigmoid')(x)
        model = tf.keras.Model(inputs, outputs)
        model.summary()
        return model

    def train(self, X_data, epochs, batch_size=128):
        X_data = self.MinMaxNorm(tf.math.log1p(X_data))
        valid = np.ones((batch_size), 1)
        fake = np.zeros((batch_size), 1)
        for epoch in range(epochs):
            '''
            Train Discriminator
            '''
            # Select random half of X_data
            idx = np.random.randint(0, X_data.shape[0], batch_size)
            bulk = X_data[idx]
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
            print ("%d [D loss: %f, acc.: %.2f%%] [G loss: %f]" % (epoch, d_loss[0], 100*d_loss[1], s_loss))
