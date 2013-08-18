import templates
import process

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# general info
detector = 'H1'
psr = 'J0534+2200'


# get data
data = process.Data(detector, psr)
data.get()

finehet = data.finehet
t = data.finehet.index

fig = plt.figure()

# generate noise
# for nlevel in [1e-30, 1e-29, 1e-28, 1e-27, 1e-26, 1e-25, 1e-24, 1e-23, 1e-22]:
#     print nlevel
d = np.random.normal(scale=nlevel, size=len(t)) + 1j*np.random.normal(scale=nlevel, size=len(t))

d = finehet
signal = templates.Signal(detector, psr, 'p', t)

# for h in [1e-27, 1e-26, 1e-25, 1e-24, 1e-23, 1e-22, 1e-21, 1e-20, 1e-19]:
for d_inj in np.linspace(1.-.001, 1.001, 10):
    
    # injection
    h = 1e-23
#     d_inj = 1.
    sweep = .001


    psi  = signal.response.src.param['POL']
    iota = signal.response.src.param['INC']

    d += h * signal.simulate(d_inj, psi, iota, phase=0.0)

    # search
    dsrch = np.linspace(1.-sweep, 1.+sweep, 101)

    res = pd.Series(index=dsrch)

    for delta in dsrch:

        s = signal.simulate(delta, psi, iota, phase=0.0)  

        res[delta] = abs(np.vdot(s, d))

        # plot
    res.plot()
    
    plt.xlim(min(dsrch), max(dsrch))
    plt.title('$\delta_{inj}=' + str(d_inj) + '$\n($h_0=' + str(h) + '$ LIGO noise)')
    plt.xlabel('$\delta_{rec}$')
    plt.ylabel('Signal - data cross correlation')
    plt.savefig('files/plots/distribution/actual_noise/varying_d/d' + str(d_inj) +'.png', bbox_inches='tight')
    plt.close()