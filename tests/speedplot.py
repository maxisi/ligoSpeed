import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sys import argv

_, kind = argv

# load results
try:
    f = pd.HDFStore('files/analysis/results/speed', 'r')
    
    results = f['crab']
    
finally:
    f.close()

# setup axes
h0 = [1e-27, 1.112e-24, 2.223e-24, 3.334e-24, 4.445e-24, 5.556e-24, 6.667e-24, 7.778e-24, 8.889e-24, 1e-23]
results.index = h0

dinj = results.columns.tolist()


if kind  == '2D':

    # plot 2D
    plt.figure()

    for h in h0:
        plt.plot(dinj, results.ix[h], '+', label=str(h))
    
    plt.legend(numpoints=1,ncol=2)

    plt.xlim(min(dinj), max(dinj))
    plt.ylim(min(dinj), max(dinj))

    plt.ylabel('Recovered $\delta$')
    plt.xlabel('Injected $\delta$')

    plt.title('Recovered $\delta=c/c_{gw}$ for different injection strengths')
    
    plt.show()

#     plt.savefig('files/plots/detection/crab_H1_with_noise_2D_2.png', bbox_inches='tight')

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
    
    plt.savefig('files/plots/detection/crab_H1_with_noise_3D_2.png')
        
else:
    print 'Error: supply kind 2D or 3D.'