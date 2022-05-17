import os
import time
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow.keras import backend as K
import numpy as np


'''
References:
github.com/eriklindernoren/Keras-GAN
https://www.tensorflow.org/tutorials/quickstart/advanced
https://www.tensorflow.org/tutorials/generative/dcgan
'''


BCE_loss = tf.keras.losses.BinaryCrossentropy(from_logits=True, name='train_loss')
Optmzr_sim = tf.keras.optimizers.Adam(0.0002, 0.5)
Optmzr_dis = tf.keras.optimizers.Adam(0.0002, 0.5)


def MinMaxNorm(x):
    x_scaled = tf.math.divide_no_nan(
        (x - tf.math.reduce_min(x)),
        (tf.math.reduce_max(x) - tf.math.reduce_min(x)))
    return x_scaled

def discriminator_loss(real, fake):
    real_loss = BCE_loss(tf.ones_like(real), real)
    fake_loss = BCE_loss(tf.zeros_like(fake), fake)
    net_loss = real_loss + fake_loss
    return net_loss

def simulator_loss(fake):
    return BCE_loss(tf.ones_like(fake), fake)


@tf.function
def train_step(real):
    noise = tf.random.normal([1,100])

    with tf.GradientTape() as tape_sim, tf.GradientTape() as tape_dis:
        fake = simulator(noise, training=True)

        real_output = discriminator(real, training=True)
        fake_output = discriminator(fake, training=True)

        loss_sim = simulator_loss(fake_output)
        loss_dis = discriminator_loss(real_output, fake_output)

        acc_dis = accuracy(real_output, fake_output)

    grad_sim = tape_sim.gradient(loss_sim, simulator.trainable_variables)
    grad_dis = tape_dis.gradient(loss_dis, discriminator.trainable_variables)

    Optmzr_sim.apply_gradients(zip(grad_sim, simulator.trainable_variables))
    Optmzr_dis.apply_gradients(zip(grad_dis, discriminator.trainable_variables))
    print ("Step: [D loss: %f] [G loss: %f]" % (loss_dis, loss_sim))


def train(dataset, epochs):
    for epoch in range(epochs):
        start = time.time()

        for image in dataset:
            train_step(image)

        if (epoch+1) % 15 == 0:
            checkpoint.save(file_prefix = checkpoint_prefix)

        print(f'Time for epoch {epoch+1} is {time.time()-start} sec')


class DirichletSimulator(tf.keras.layers.Layer):

    def __init__(self, total_count=500):
        super(DirichletSimulator, self).__init__()
        self.total_count=total_count

    def call(self, inpt):
        alpha = inpt#.numpy().tolist() # provide as list or tfp function will run forever
        dist = tfp.distributions.DirichletMultinomial(total_count=[self.total_count]*len(alpha), concentration=alpha)
        outpt = dist.sample() # no argument to ".sample()" produces output of shape (batchsize, # of classes)
        return tf.cast(outpt, dtype=tf.int32) # "tf.slice" takes in tensor of integer values


class Subsampling(tf.keras.layers.Layer):

    def __init__(self, scdata):
        super(Subsampling, self).__init__()
        # initialize one layer for each cell type
        self.scdata=scdata

    def call(self, inpt):
        subset = tf.slice(tf.random.shuffle(self.scdata), begin=[0,0], size=[inpt.numpy(), -1])
        outpt = tf.reduce_sum(subset, axis=0)
        return outpt


class Build_Simulator(tf.keras.Model):

    def __init__(self, scdata, n_cellspersample=500):
        super(Build_Simulator, self).__init__()
        self.scdata = scdata
        self.n_celltypes = len(scdata)
        self.Dense = tf.keras.layers.Dense(10, activation='relu')
        self.Dirichlet = DirichletSimulator(total_count=n_cellspersample)
        self.Add = tf.keras.layers.Add()
        self.Log1p = tf.keras.layers.Lambda(lambda t: tf.math.log1p(t))
        self.Norm = tf.keras.layers.LayerNormalization(center=True, scale=True)
        self.MinMax = tf.keras.layers.Lambda(lambda t: MinMaxNorm(t))

    def call(self,x):
        x=self.Dense(x)
        x=self.Dirichlet(x)
        x_=[]
        for c in range(self.n_celltypes):
            x_.append(Subsampling(self.scdata[c])(x[0][c]))
        x=self.Add(x_)
        x=self.Log1p(x)
        x=self.Norm(x)
        x=self.MinMax(x)
        return x


class Build_Discriminator(tf.keras.Model):

    def __init__(self):
        super(Build_Discriminator, self).__init__()
        self.Dense1 = tf.keras.layers.Dense(64)
        self.Activation = tf.keras.layers.LeakyReLU(alpha=0.2)
        self.Dense2 = tf.keras.layers.Dense(32)
        self.Output = tf.keras.layers.Dense(1, activation='sigmoid')

    def call(self,x):
        x=self.Dense1(x)
        x=self.Activation(x)
        x=self.Dense2(x)
        x=self.Activation(x)
        x=self.Output(x)
        return x
