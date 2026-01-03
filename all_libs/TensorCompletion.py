# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 20:55:49 2023

@author: SmartDATALab
"""

import torch
import tntorch as tn
import numpy as np
import scipy.io as sci
import matplotlib.pyplot as plt

data = sci.loadmat('wave_data_incomplete.mat')['wave_data']



im = data
#plt.imshow(im, cmap='gray', vmin=im.min(), vmax=im.max())
#plt.show()

#P = im.shape[0]*im.shape[1]
#Q = int(P/10)
#print('We will keep {} out of {} pixels'.format(Q, P))
#X = np.unravel_index(np.random.choice(P, Q), im.shape)  # Coordinates of surviving pixels
y = torch.Tensor(im)  # Grayscale values of surviving pixels




t = tn.rand(im.shape, ranks_tt=20, requires_grad=True)

def loss(t):
    return tn.relative_error(y, t), tn.normsq(tn.partialset(t, order=2))*1e-4
tn.optimize(t, loss)




plt.imshow(t[:,:,1].numpy())
plt.show()