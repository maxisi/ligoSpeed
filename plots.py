import templates
import matplotlib.pyplot as plt
import sidereal as sd
import paths

import numpy as np
import pandas as pd
import os


# speed
def roemer(days, step=60, detector='H1', psr='J0534+2200'):
    t = np.array([x+ 630720013 for x in range(0, int(days*sd.ss), step)])
    s = templates.Signal(detector, psr, 'p', t)
    
    plt.figure()
    s.roemer.plot()
    
    plt.title('Roemer delay between ' + detector + ' and ' + psr + ' over ' + str(days) + ' sidereal days')
    plt.xlabel('GPS time (s)')
    plt.ylabel('$\Delta_R$ (s)')
    
    plt.xlim(min(s.roemer.index), max(s.roemer.index))
    
    p = paths.plots + 'roemer/'

    try:
        os.makedirs(p)
    except:
        pass
    
    plt.savefig(p + 'roemer_' + detector + '_' + psr + '_' + str(days) + 'days', bbox_inches='tight')
    plt.close()