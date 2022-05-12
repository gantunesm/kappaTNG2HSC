import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool

zs_tng = array([0.034, 0.07 , 0.105, 0.142, 0.179, 0.216, 0.255, 0.294, 0.335,
       0.376, 0.418, 0.462, 0.506, 0.552, 0.599, 0.648, 0.698, 0.749,
       0.803, 0.858, 0.914, 0.973, 1.034, 1.097, 1.163, 1.231, 1.302,
       1.375, 1.452, 1.532, 1.615, 1.703, 1.794, 1.889, 1.989, 2.094,
       2.203, 2.319, 2.44 , 2.568])

zs =['z1','z2','z3','z4' ,'singlez']


ng = 1024  # number of grids
theta = 5.0  # opening angle in deg
 
# From Jia Liu and Ken Osato codes to read the data:

def kappa_gen(LP, run, iz):
    fname = "/global/cscratch1/sd/jialiu/kappaTNG/kappaTNG-Hydro/LP%03d/run%03d/kappa%02d.dat"%(LP, run, iz)
       #kappaTNG-Hydro
    with open(fname, 'rb') as f:
        dummy = np.fromfile(f, dtype="int32", count=1)
        kappa = np.fromfile(f, dtype="float", count=ng*ng)
        dummy = np.fromfile(f, dtype="int32", count=1)
    return kappa.reshape((ng, ng))


def cosmos_kappa_gen(LP, run, weights):
    kappa_cosmos = np.zeros(shape=(ng, ng))
    for iz in range(1,41): ##### redshift counter
        kappa_cosmos += weights[iz-1]*kappa_gen(LP, run, iz)
    return kappa_cosmos  


def mass_product (counter):
       '''counter=(LP-1)*100+run'''
       run = counter%100
       if run ==0:
              run = 100
       LP = int((counter - run)/100)+1
       print(LP,'LP')
       print(run,'run')
       for zz in range(len(zs)):
              fn = '/global/cscratch1/sd/gmarques/tng/kappaTNG-Hydro_HSCzs/'+zs[zz]+'/kappa%i.npy' % (counter)
              ww = np.load('/global/homes/g/gmarques/hsc_ng/for_baryons/weights_pz/pz_'+zs[zz]+'_zstng_normed.npy')
              kappa = cosmos_kappa_gen(LP, run, ww)
              std_kappa = kappa.std()
              fn2 = '/global/cscratch1/sd/gmarques/tng/kappaTNG-Hydro_HSCzs/'+zs[zz]+'/std_kappa%i.npy' % (counter)
              np.save(fn2, std_kappa)
              np.save(fn, kappa)




pool = Pool() 
out = pool.map(mass_product, np.arange(2000, 8000))
pool.close()


 