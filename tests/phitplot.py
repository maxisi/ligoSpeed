import templates
import matplotlib.pyplot as plt
import numpy as np
import sidereal as sd

reload(templates)

days=1
t = np.array([x+ 630720013 for x in range(0, int(days*sd.ss), 60)])

plt.figure()

det = 'H1'

kind = 'GR'
pdif = 'p'

phi0 = 5.

for d in np.linspace(.999,1.001, 5):

    sig = templates.Signal(det, 'J0534+2200', pdif, t)

    s = sig.phase_ev(d, phi0)
    
    s.plot(label=str(d))
    
plt.legend(numpoints=1)
    
plt.title('Crab ' + det + 'phase evolution for ' + kind + pdif + ' varying $c/c_{gw}$ and $\phi_0=$' + str(phi0))

plt.show()