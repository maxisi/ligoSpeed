import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sys import argv

_, kind = argv

# load results
path = 'files/analysis/results/small_delta_gauss'
try:
    f = pd.HDFStore(path, 'r')
    results = f['crab']
    
finally:
    f.close()

# results = pd.read_pickle('files/analysis/results/highSignal')

# setup axes
h0 = results.index

dinj = results.columns.tolist()


if kind  == '2D':

    # plot 2D
    plt.figure()

    for h in h0[:7]:
        plt.plot(dinj, results.ix[h], label=str(h))
    

    plt.xlim(min(dinj), max(dinj))
#     plt.ylim(1.-.016, 1.011)#.98, 1.+.101)
#     
#     # search limits
#     plt.plot(dinj, [1.+.1]*len(dinj), '--', color='.7')
#     plt.plot(dinj, [1.-.015]*len(dinj), '--', color='.7')

    plt.legend(numpoints=1,ncol=2, loc=4)

    plt.ylabel('Recovered $\delta$')
    plt.xlabel('Injected $\delta$')

    plt.title('Recovered $\delta=c/c_{gw}$ for different injection strengths\n(S6 H1 Crab noise std ~1e-23)')
    
    plt.show()

#     plt.savefig('files/plots/detection2/crab_H1_LIGO_noise_h-25-18_2D-line_smalld.png', bbox_inches='tight')

elif kind == '3D':
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for h in h0:
        x = [h] * len(dinj)
        y = dinj
        z = results.ix[h].tolist()

        ax.scatter(x, y, z)

        xLabel = ax.set_xlabel('$h_0$')
        yLabel = ax.set_ylabel('Injected $\delta$')
        zLabel = ax.set_zlabel('Recovered $\delta$')
        

    plt.title('Recovered $\delta=c/c_{gw}$')
    
#     plt.savefig('files/plots/detection2/crab_H1_with_noise_3D_2.png')
    
    plt.show()

        
else:
    print 'Error: supply kind 2D or 3D.'