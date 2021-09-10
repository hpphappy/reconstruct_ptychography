import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('pdf')
import h5py
import os
from util import *
import dxchange
import numpy as np

path = 'cell/NMF_ptycho_recon'

save_path = os.path.join(path, 'comparison')
if not os.path.exists(save_path):
    os.makedirs(save_path)

grid_delta = np.load('cell/phantom/grid_delta.npy')
print('dimension of the sample = ' +', '.join(map(str,grid_delta.shape)))
grid_delta = np.squeeze(grid_delta)

n_sample_pixel = grid_delta.shape[0]*grid_delta.shape[1]

energy_kev = 5
k = (2*np.pi)/(1.24/energy_kev) #k in the unit of nm^-1
print(k)

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family': 'serif', 'serif': ['Helvetica'], 'weight': 'normal', 'size': 12}
plt.rc('font', **fontProperties)

n_s_ls = [1, 10, 30, 50, 100, 300, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]

# spec = gridspec.GridSpec(1, 2 , width_ratios=[7, 0])
smse_normal_ls = []
fig, ax = plt.subplots(figsize=(8,5))
for i, n_s in enumerate(n_s_ls):
    obj_dir = os.path.join(path, 'n2e7_nei' + str(n_s))
    obj = dxchange.read_tiff(os.path.join(obj_dir, 'delta_ds_1.tiff'))
    obj = obj[:,:,0]
    smse_normal = (k**2)*np.sum((grid_delta - obj)**2)/n_sample_pixel
    smse_normal_ls.append(smse_normal)

np.savez(os.path.join(save_path, "mse"), n_s_ls, smse_normal_ls)
ax.plot(n_s_ls, smse_normal_ls, '-s', color='#C000FF', linewidth=2, markerfacecolor='none', markeredgecolor='#C000FF', markeredgewidth=2, alpha=0.5, label = 'PCA')
ax.set_xticks(np.arange(0,4000,500))


Fontsize = 12
ax.set_xlabel('S', fontsize=Fontsize)
ax.set_ylabel('MSE', fontsize=Fontsize)
ax.legend(loc=1,fontsize=Fontsize,ncol=1)
# ax.set_xscale('log')
ax.set_yscale('log')

plt.savefig(os.path.join(save_path, 'MSE.pdf'), format='pdf')
#plt.savefig(os.path.join(save_path, 'SMSE.png'))
fig.clear()
plt.close(fig)