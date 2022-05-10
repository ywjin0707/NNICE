# -*- coding: utf-8 -*-
"""
Created on Thu May  5 19:58:48 2022

@author: ywjin0707
"""

import tensorflow as tf


def pearson_cor(targets, logits):
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
    r = tf.math.divide_no_nan(r_num, r_den)
    r = tf.maximum(tf.minimum(r, 1.0), -1.0)
    return r

def tilted_loss(targets, logits, q, alpha=0.9):
    r = pearson_cor(targets, logits)
    rmse = tf.math.sqrt(tf.math.reduce_mean(tf.math.squared_difference(targets, logits)))
    e = (alpha*(1 - r)) + (1 - alpha)*rmse
    return tf.math.maximum(q * e, (q - 1)*e)

