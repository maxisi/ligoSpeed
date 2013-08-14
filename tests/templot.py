import templates
import matplotlib.pyplot as plt
import numpy as np
import sidereal as sd

import cmath

reload(templates)

days=1
t = np.array([x+ 630720013 for x in range(0, int(days*sd.ss), 60)])


det = 'H1'

kind = 'GR'
pdif = 'm'

phi0 = np.pi/4.

dd = .5
drange = np.linspace(1-dd,1+dd, 5)

plt.figure()

# for d in drange:
# 
#     sig = templates.Signal(det, 'J0534+2200', pdif, t)
# 
#     s = sig.simulate(d, sig.response.src.param['POL'], sig.response.src.param['INC'], phase=phi0)
#     
#     abs(s).plot(label=str(d))
#     
# plt.legend(numpoints=1)
#     
# plt.title('Crab ' + det + ' response amplitude for ' + kind + pdif + ' varying $c/c_{gw}$ and $\phi_0=$' + str(phi0))
# plt.xlabel('GPS time (s)')
# plt.ylabel(' Response amplitiude')
# plt.savefig('files/plots/response/p0'+ str(phi0) + '_crab_' + kind + pdif + '_' + det + '_response_A_'+ str(phi0) +'.png', bbox='tight')
# 
# 
# plt.figure()
# 
# for d in drange:
# 
#     sig = templates.Signal(det, 'J0534+2200', pdif, t)
# 
#     s = sig.simulate(d, sig.response.src.param['POL'], sig.response.src.param['INC'], phase=phi0)
#     
#     s.imag.plot(label=str(d))
#     
# plt.legend(numpoints=1)
#     
# plt.title('Crab ' + det + ' response (Im) for ' + kind + pdif + ' varying $c/c_{gw}$ and $\phi_0=$' + str(phi0))
# 
# plt.xlabel('GPS time (s)')
# plt.ylabel(' Response (Re)')
# plt.savefig('files/plots/response/p0'+ str(phi0) + '_crab_' + kind + pdif + '_' + det + '_response_Im_'+'.png', bbox='tight')
# 
# 
# 
# plt.figure()
# 
# for d in drange:
# 
#     sig = templates.Signal(det, 'J0534+2200', pdif, t)
# 
#     s = sig.simulate(d, sig.response.src.param['POL'], sig.response.src.param['INC'], phase=phi0)
#     
#     s.real.plot(label=str(d))
#     
# plt.legend(numpoints=1)
#     
# plt.title('Crab ' + det + ' response (Re) for ' + kind + pdif + ' varying $c/c_{gw}$ and $\phi_0=$' + str(phi0))
# 
# plt.xlabel('GPS time (s)')
# plt.ylabel(' Response (Re)')
# plt.savefig('files/plots/response/p0'+ str(phi0) + '_crab_' + kind + pdif + '_' + det + '_response_Re_'+ str(phi0) +'.png', bbox='tight')


plt.figure()

for d in drange:

    sig = templates.Signal(det, 'J0534+2200', pdif, t)

    s = sig.simulate(d, sig.response.src.param['POL'], sig.response.src.param['INC'], phase=phi0)
    
    s.map(cmath.phase).plot(label=str(d))
    
plt.legend(numpoints=1, loc=4)
    
plt.title('Crab ' + det + 'response amplitude for ' + kind + pdif + ' varying $c/c_{gw}$ and $\phi_0=$' + str(phi0))

plt.xlabel('GPS time (s)')
plt.ylabel(' Response (phase)')
plt.savefig('files/plots/response/p0'+ str(phi0) + '_crab_' + kind + pdif + '_' + det + '_response_P_'+ str(phi0) +'.png', bbox='tight')



plt.close()