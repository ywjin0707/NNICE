# -*- coding: utf-8 -*-
"""
Created on Thu May  5 20:51:41 2022

@author: ywjin0707
"""
sc.pp.normalize_total(mydata, target_sum = 1e4)
sc.pp.log1p(mydata)
sc.pp.highly_variable_genes(mydata)
mydata = mydata[:, mydata.var.highly_variable]


def simulate_bulk(datasets: list, n: int = 10000, c: int = 500, sf: int = 100):
    k = len(datasets)
    Nprop = np.random.dirichlet(alpha=[1]*k, size=n)
    for prop in Nprop:
        bulk_sample = [ds[ds.obs.index.isin(np.random.randint(low=0, high=ds.n_obs, size=prop).tolist())] for ds in datasets]
        bulk_sample = ad.concat(bulk_sample)
        