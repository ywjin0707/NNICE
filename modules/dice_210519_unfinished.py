# -*- coding: utf-8 -*-
"""
Created on Mon May  3 23:11:04 2021

@author: yw_ji
"""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = 2
import datetime
import tensorflow as tf
    
    
class DICE(object):
    def __init__(
            self,
            model_dir,
            model_name,
            batch_size=100,
            learn_rate=1e-3,
            epochs=1000,
            lmbda=0.0001,
            alpha=0.8
    ):
        self.model_dir = model_dir
        self.model_name = model_name
        self.batch_size = batch_size
        self.learn_rate = learn_rate
        self.epochs = epochs
        self.lmbda = lmbda
        self.alpha = alpha
        
        self.datetime = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        self.quantiles = [0.1,0.25,0.5,0.75,0.9]
        
        self.train_data = None
        self.valid_data = None
        self.celltypes = None
        self.n_celltypes = None
        self.num_sig_genes = None
    
    def build_model(self):
        model = tf.keras.Sequential()
        model.add(tf.keras.layers.Dense(1, activation='linear', bias_constraint=tf.keras.constraints.NonNeg(),
                                        kernel_regularizer=tf.keras.regularizers.l2(self.lmbda)))
        return model
    
    def pearson_cor(self, targets, logits):
        mx = tf.reduce_mean(logits)
        my = tf.reduce_mean(targets)
        xm, ym = logits - mx, targets - my
        r_num = tf.reduce_sum(input_tensor=tf.multiply(xm, ym))
        r_den = tf.sqrt(
            tf.multiply(
                tf.reduce_sum(input_tensor=tf.square(xm)),
                tf.reduce_sum(input_tensor=tf.square(ym)),
            )
        )
        r = tf.divide(r_num, r_den)
        r = tf.maximum(tf.minimum(r, 1.0), -1.0)
    return r
    
    def loss_obj(self, targets, logits, q):
        r = self.pearson_cor(targets, logits)
        rmse = tf.math.sqrt(tf.math.reduce_mean(tf.math.squared_difference(targets, logits)))
        e = (self.alpha*(1 - r)) + (1 - self.alpha)*rmse
        return tf.math.reduce_max(q * e, (q - 1)*e)

    def load_data(self, input_path):
        x_data = pd.read_csv(self.input_path, index_col=0)
        y_data = pd.read_csv(self.input_path, index_col=0)
        
    def compile_fit_save(self, model, c, q):
        model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=self.learn_rate),
                      loss=lambda y,f: self.loss_obj(y,f,q),
                      metrics=[self.pearson_cor, tf.keras.metrics.RootMeanSquaredError()])
        log_name = f'{self.model_dir}{self.datetime}_C{c}_Q{q*100}'
        model.fit(self.train_data, 
                  batch_size=self.batch_size, 
                  epochs=self.epochs,
                  validation_data=self.valid_data,
                  callbacks=[tf.keras.callbacks.CSVLogger(filename=log_name+'.log')],
                  verbose=0)
        model.save(log_name+'.')
        res = model.evaluate(valid_data)
        print(f'Validation metrics --- Celltype:{c}  Quantile:{q*100}  Loss:{res[0]}  R:{res[1]}  RMSE:{res[2]}')
        return model
    
    def train(self,):
        for C in self.celltypes:
            pred = []
            true = []
            for Q in self.quantiles:
                


import pandas as pd
pd.read_table('C:/Users/yw_ji/Documents/MSc Thesis/DATA/scaden_input/')

import anndata as ad
x = ad.read_h5ad('C:/Users/yw_ji/Documents/MSc Thesis/DATA/scaden_input/simulated_data/data.h5ad')
x_data = x.X.astype(np.float32)
fractions = [x.obs[ctype] for ctype in x.uns["cell_types"]]
y_data = np.array(fractions, dtype=np.float32).transpose()
data = tf.data.Dataset.from_tensor_slices((self.x_data, self.y_data))
data = data.shuffle(100).repeat().batch(32)

model = tf.keras.Sequential()
model.add(SignatureMatrix(526))

model.compile(optimizer='Adam', loss='mse', metrics=[pearson_cor])

model.fit(data, steps_per_epoch=1000, epochs=10)