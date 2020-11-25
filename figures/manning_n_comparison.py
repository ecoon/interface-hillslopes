import numpy as np
from matplotlib import pyplot as plt
import colors # $ATS_SRC_DIR/tools/utils

cm = colors.cm_mapper(0,6)

def plot(mann, ax, color='b', marker=''):
    dat = np.loadtxt(f'../model_03-manning_n/manning_n-{mann}/fluxes.dat')
    dtotal = dat[0:365,1].sum()
    dmax = dat[0:365,1].max()
    ax.plot(dat[:,0], dat[:,1]/55000, '-'+marker, color=color, label=f'mann_n = {mann}')
    print(f'manning n = {mann}: discharge total = {dtotal}; max = {dmax}')

fig = plt.figure()
ax = fig.add_subplot(111)
for i, man in enumerate([0.1, 0.2, 0.5, 1, 2, 5, 10]):
    plot(man, ax, color=cm(i))

ax.set_xlabel('time [doy]')
ax.set_ylabel('discharge [m^3/day]')
ax.legend()
plt.savefig('manning_n_comparison')
plt.show()
