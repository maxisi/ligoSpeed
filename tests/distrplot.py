import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sys import argv

_, kind = argv

# load results
try:
    f = pd.HDFStore('files/analysis/results/distribution', 'r')
    
    results = f['crab']
    
finally:
    f.close()

# setup axes
drec = results.index.tolist()
dinj = results.columns


if kind  == '2D':

    # plot 2D
    plt.figure()

    for d in dinj:
        results[d].plot()
#         plt.plot(drec, results[d], '+', label=str(d))
    
#     plt.legend(numpoints=1,ncol=2)

    plt.xlim(min(drec), max(drec))
#     plt.ylim(min(drec), max(drec))

    plt.ylabel('Cross correlation')
    plt.xlabel('Recovered $\delta$')

#     plt.title(' $\delta=c/c_{gw}$ for different injection strengths')
    
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