from sys import exit
import matplotlib.pyplot as plt
import os
import datetime
from time import time

import numpy as np
import pandas as pd
import random
import math
import cmath
from scipy import signal
import random

import templates
import sidereal as sd
import paths
import psrplot

reload(sd)
reload(templates)
reload(psrplot)


def het(vector, f, *arg):
    '''
    Heterodynes vector at frequencies f. Preferred input for vector is series indexed
    over time; f can be a list or an array. Returns DataFrame.
    '''
    print 'Ready to heterodyne.',
    
    if len(arg)==0:
        try:
            t = vector.index.tolist()
        except AttributeError:
            print 'ERROR: no time vector for heterodyne.'
            exit(0)
            
    elif len(arg)==1:
        t = arg[0]
        
    else:
        print 'ERROR: het needs input time or indexed vector, not %d extra arguments.' % len(arg)
        exit(0)
    
    temp = np.exp(2*np.pi*1j*np.multiply.outer(f, t))
    
    try:
        template = pd.DataFrame(temp, columns=t)
    except ValueError:
        template = pd.Series(temp, index=t)
    
    rh = vector*template
    print 'Done'
    return rh.T


class Data(object):
    '''
    Holds original data and related information.
    '''
    def __init__(self, detector, psr):        
        self.detector = detector
        self.det = sd.detnames(detector)
        self.psr = psr
        
        # data info
        self.datadir = paths.importedData + self.psr + '_' + self.detector + '.hdf5'
        self.seedname = 'finehet_' + self.psr + '_' + self.detector


    def imp(self):
        '''
        Return DF with original data (col: PSR; index: t). Assuming execution on ATLAS.
        '''
        
        struct = '/data' + self.detector + '/' + self.seedname
        pathOptions = [
                     paths.originalData + struct,
                     paths.originalData + '/' + self.psr + '_' + self.detector + struct
                     ]
        
        try:
            d = pd.HDFStore(self.datadir, 'w')
                
            for p in pathOptions:
                try:
                    dParts = pd.read_table(p, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                except IOError:
                    pass
            
            # check file was found
            try:
                dParts
            except NameError:
                print 'Could not find %s data for PSR %s in:\n' % (self.detector, self.psr),
                for p in pathOptions:
                    print '\t%(p)s' % locals()
                    
                print 'Should I...'
                print '\t1. Provide path\n\t2. Abort'
                opt = raw_input('>')
            
                if opt==1:
            
                    try:
                        psrpath = raw_input('Enter path: ')
                        dParts = pd.read_table(psrpath, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                
                    except IOError:
                        print "File not found. One more try?"
                        psrpath = raw_input('Enter path: ')
                        dParts = pd.read_table(psrpath, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
            
                else:
                    print 'Could not find %(detector)s data for PSR %(psr)s.' % locals()
                    print 'Exiting at analysis/process ln 77'
                    exit()

            self.finehet = dParts['Re']+dParts['Im']*1j
            d[self.psr] = self.finehet
            
        finally:
            d.close()


    def get(self):
        '''
        Retrieves original heterodyned data for pulsars in list.
        Imports data from M.Pitkin if necessary.
        '''
        
        try:
            d = pd.HDFStore(self.datadir, 'r')
        
            try:
                self.finehet = d[self.psr]
            except KeyError:
                # file is empty or is corrupted
                d.close()
                self.imp()
            else:
                d.close()
        
        except IOError:
            self.imp()
                

class Background(object):
    '''
    Manages background files for a given detector and source: gets and creates.
    '''
    def __init__(self, detector, psr, freq, filesize=100):
        # data
        self.seed = Data(detector, psr)
        self.seed.get()
        
        # background info
        self.freq = freq
        self.filesize = filesize      # number of series per file. Adjust!
        
        self.nsets = int(len(freq)/filesize)
        if self.nsets<1:
            self.nsets = 1
            
        self.fset = {n : freq[n*filesize:min(len(freq),(n+1)*filesize)] for n in range(self.nsets)}
        
        # storing info
        self.dir = paths.rhB + self.seed.det + '/' + psr + '/'
        self.name = 'back_' + psr + '_' + self.seed.det + '_'
        self.path = self.dir + self.name
    
    def writelog(self):
        now = datetime.datetime.now()
        comments = '# ' + self.seed.detector + '\n# ' + self.seed.psr + '\n# ' + str(now) + '\n'
        fileinfo = 'nsets\tfilesize\n' + str(self.nsets) + '\t' + str(self.filesize)
        
        try:
            f = open(self.dir + 'log.txt', 'w')
            f.write(comments + fileinfo)
        finally:
            f.close()  
            
               
    def create(self):
        '''
        Re heterodynes and saves data at frequencies f. Number of heterodynes is determined by
        f and data can be for more than one pulsar
        '''
        
        # create background directory
        try:
            os.makedirs(self.dir)
        except OSError:
            pass
            
        # create background files
        for n in range(self.nsets):
            path = self.dir + self.name + str(n)
            
            try:
                rh = pd.HDFStore(path, 'w')
                rh[self.seed.psr] = het(self.seed.finehet, self.fset[n])
            finally:
                rh.close()
        
        self.writelog()

    def get(self):
        '''
        Checks background required for search exits and creates it if needed.
        Returns filename list.
        '''
        # read log
        try:
            readme = pd.read_table(self.dir + 'log.txt', sep='\s+', skiprows=3)
            log_nfiles = readme['nsets'].ix[0]
            log_filesize = readme['filesize'].ix[0]
            log_nfreq = log_nfiles * log_filesize
            
            # get actual number of background files in directory
            files = [name for name in os.listdir(self.dir) if 'back' in name]
            nfiles = len(files)
            
            if nfiles!=log_nfiles or log_nfreq!=len(self.freq) or log_filesize!=self.filesize:
                self.create()
        except IOError:
            # no log found
            self.create()


class Results(object):
    '''
    Holds search results and contains methods to save them.
    '''
    def __init__(self, detector, psr, dinj=[], hinj=[], pdif='p', extra_name=''):
        # system
        self.detector = detector
        self.psr = psr
        
        # search
        self.methods = ['GR']
        
        # injection
        self.dinj = dinj
        self.hinj = hinj
        
        self.kind = 'GR'
        
        self.pdif = pdif
        
        # containers
        self.delta =  pd.DataFrame(columns = dinj, index = hinj)
                
        # saving
        self.extra_name = extra_name
        self.dir = paths.results + self.detector + '/' + self.psr + '/' 
        self.name = 'drec_' + self.psr + '_' + self.detector + '_' + sd.phase2(pdif) + '_' + extra_name
        self.path = self.dir + self.name
        
        self.issaved =  False
        
    def save(self):
        
        try:
            os.makedirs(self.dir)
        except OSError:
            pass
        
        try:
            os.remove(self.path)
        except:
            pass
                
        try:
            self.delta.to_pickle(self.path)
            print '\nResults saved to:'
            print self.path
            self.issaved = True
        except IOError:
            print 'Failed to save!'
            print self.path 
            
        
    def load(self):
        self.delta = pd.read_pickle(self.path)


    def plot(self, kind='2D', extra_name='', save=True, style='-', lgloc=4):
        
        header = 'Recovered $\delta=c/c_{gw}$ for different injection strengths ' + extra_name
        
        results = self.delta
        h0 = results.index
        dinj = results.columns.tolist()        

        if kind  == '2D':

            # plot 2D
            plt.figure()

            for h in h0:
                plt.plot(dinj, results.ix[h], style, label=str(h))
    

            plt.xlim(min(dinj), max(dinj))
            plt.ylim(results.min().min()-.001, results.max().max()+.001)
    
            plt.legend(numpoints=1,ncol=2, loc=lgloc)

            plt.ylabel('Recovered $\delta$')
            plt.xlabel('Injected $\delta$')


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
        
        else:
            print 'Error: supply kind 2D or 3D.'
            raise NameError
        
        plt.title(header)

        plt.show()
        
        if save:
            pltdir = paths.plots + self.detector + '/detection/'
            pltname = 'delta_' + self.detector + '_' + self.psr + '_' + sd.phase2(self.pdif) + '_' + kind + '_' + extra_name
            save_to = pltdir + pltname
            print 'Plot saved to:\n %(save_to)s' % locals()
        
            try:
                os.makedirs(pltdir)
            except OSError:
                pass
            
            plt.savefig(save_to, bbox_inches='tight')
            
            plt.close()
        


class Delta(object):
    
    def __init__(self, detector='H1', psr='J0534+2200', ninj=8, h0=-25, pdif='p', nd=20, sweep=101, dinj0=1., dinjd=.01):
        # system info
        self.detector = detector
        self.psr = psr
                
        # data info
        self.data = Data(detector, psr)
        self.data.get()
        
#         frange = [1.0e-7, 1.0e-5]
#         self.freq = np.linspace(frange[0], frange[1], nfreq)
#         print 'Getting background.'
#         self.background = Background(detector, psr, self.freq, filesize=100)
#         self.background.get()
#         
        self.t = self.data.finehet.index
        
        sigma = Sigma(self.detector, self.psr, self.data.finehet)
        self.sg = sigma.std
        
        # injection strengths
        self.hinj = [10**(p) for p in range(h0, h0+ninj)]
        
        # injection deltas
        self.dinj = np.linspace(dinj0, dinj0+dinjd, nd)
        
        # injection kind   
        self.pdif = pdif
        self.injkind = 'GR'
        
        # set up injection
        self.signal = templates.Signal(detector, psr, pdif, self.t)
        self.setranges([])
        
        # search
        self.dsrch = np.linspace(dinj0-.001, dinj0+dinjd+.001, sweep)
        
                
    def setranges(self, rangeparam):
        src = self.signal.response.src
        
        if rangeparam ==[]:
            print 'All parameters fixed.'
        else:
            print 'Ranging over: ' + str(rangeparam)
        
        
        if 'psi' in rangeparam or rangeparam=='all':
            self.pol_range = [
                            src.param['POL'] - src.param['POL error'],
                            src.param['POL'] + src.param['POL error']
                            ]
        else:
            self.pol_range = [src.param['POL'], src.param['POL']]
        
        if 'iota' in rangeparam or rangeparam=='all':   
            self.inc_range = [
                            src.param['INC'] - src.param['INC error'],
                            src.param['INC'] + src.param['INC error']
                            ]
        else:
            self.inc_range = [src.param['INC'], src.param['INC']]

        if 'phi0' in rangeparam or rangeparam=='all':                     
            self.phi0_range = [0., np.pi/2.]
        else:
            self.phi0_range = [0., 0.]
    
    
    def recover(self, noisetype='s6', extra_name=''):
        
        # messages
        print 'Recovering delta in Gaussian noise.'
        print '\nDelta injections:'
        print self.dinj
        print '\nDelta search:'
        print self.dsrch
        print '\nInjection strengths:'
        print self.hinj  

        # results
        self.results = Results(self.detector, self.psr, self.dinj, self.hinj, pdif = self.pdif, extra_name=noisetype + extra_name)
            
        # select injection psi and iota
        psi  = random.uniform(self.pol_range[0], self.pol_range[1])
        iota = random.uniform(self.inc_range[0], self.inc_range[1])  
        phi0 = 0.0
        
        # form noise
        if noisetype=='gaussian':
            nlevel = 1.2e-23
            noise = np.random.normal(scale=nlevel, size=len(self.t)) + 1j*np.random.normal(scale=nlevel, size=len(self.t))
            print '\nUsing gaussian noise std: ' + str(nlevel)
            
        elif noisetype == 's6':
            noise = self.data.finehet
            print '\nUsing LIGO S6 noise.'  
        
        # compare sigmas
        plt.figure()
        self.sg.plot(style='+')
        plt.plot(self.t, [np.std(noise)]*len(self.t), color='r')
        plt.title('Standard deviation comparison: S6 data vs. fabricated noise')
        plt.xlabel('GPS time')
        plt.ylabel('$\sigma$')
        plt.savefig('files/plots/sigmacomp.png', bbox_inches='tight')
        plt.show()
        plt.close()
        
        print 'Creating templates.'
        
        # create injection signals
        sinj = {d_inj : self.signal.simulate(d_inj, psi, iota) for d_inj in self.dinj}

        # create search signals
        ssrch = {d_srch : self.signal.simulate(d_srch, psi, iota) for d_srch in self.dsrch}
        
        
        print '\nProcessing:'
        
        # process
        for h in self.hinj:
            print h,
        
            for d_inj in self.dinj:
            
                # inject signal
                d = h * sinj[d_inj] + noise
                
                # search
                cor = []
                for d_srch in self.dsrch:
                    cor += [abs(np.vdot(ssrch[d_srch], d))]
                
                self.results.delta[d_inj][h] = self.dsrch[np.argmax(cor)]
        
        print '\n'
               
        self.results.save()


    def confidence(self, noisetype='s6', h=1e-23, loops=100):
        
        # results
        self.results = pd.DataFrame(index=range(loops), columns=self.dinj)
        
        # select injection psi and iota
        psi  = random.uniform(self.pol_range[0], self.pol_range[1])
        iota = random.uniform(self.inc_range[0], self.inc_range[1])  
                
        print 'Obtaining uncertainty range for h0 = ' + str(h)
        self.hinj = [h] * loops

        nlevel = 1.2e-23
        
        print 'Creating templates.'
        
        # create injection signals
        sinj = {d_inj : self.signal.simulate(d_inj, psi, iota) for d_inj in self.dinj}

        # create search signals
        ssrch = {d_srch : self.signal.simulate(d_srch, psi, iota) for d_srch in self.dsrch}
        
        # analyze
        for n in range(loops):

            print n
            
            noise = np.random.normal(scale=nlevel, size=len(self.t)) + 1j*np.random.normal(scale=nlevel, size=len(self.t))
        
            for d_inj in self.dinj:
            
                # inject signal
                d = h * sinj[d_inj] + noise
                
                # search
                cor = []
                for d_srch in self.dsrch:
                    cor += [abs(np.vdot(ssrch[d_srch], d))]
                
                self.results[d_inj][n] = self.dsrch[np.argmax(cor)]
        
        # histogram
        self.variation = self.results - self.dinj
        
        plt.figure()
        
        self.variation.hist(sharex=True, sharey=True, xrot=90., figsize=(12.,8.))
        
        path = paths.plots + 'confidence/' + self.detector + '/'
        
        try:
            os.makedirs()
        except:
            pass
        
        plt.suptitle('$\delta$ recovery residuals for $h_0=' + str(h) + '$ in gaussian noise for ' + str(loops) + ' iterations\n(a residual of 0 means $\delta{inj}=\delta_{rec}$')
        
        plt.savefig(path + 'residuals_' + str(h), bbox_inches='tight')
        plt.close()
        
        plt.figure()
        
        # overall distance
        self.distance = np.sqrt((self.variation**2).sum(axis=1))
        self.distance.hist(bins=10)
        plt.title('Distance between $\delta_{rec}$ and $\delta_{inj}$ for $h_0=' + str(h) + '$\nin gaussian noise over ' + str(loops) + ' iterations and $\delta_{inj} \in [' + str(min(self.dinj)) + ',' + str(max(self.dinj)) + ']$')
        plt.xlabel('Euclidean distance')
        plt.ylabel('Count')
        plt.savefig(path + 'distance_' + str(h), bbox_inches='tight')
        plt.close()
                    

# class Distribution(InjSearch):
#     
#     def __init__(self):
#         super(Distribution, self).__init__('H1', 'J0534+2200', 10, 0., 50)
#         
#     def analyze(self):
# 
#         print 'Getting matched filter distribution.'
#     
#         # search info
#         search = templates.Signal(self.detector, self.psr, self.pdif, self.t)
# 
#         # results
#         self.results = pd.DataFrame(columns=self.dinj, index=self.hinj)
#             
#         # loop over files
#         for n in range(self.background.nsets):
#             
#             try:
#                 back_file = pd.HDFStore(self.background.path + str(n), 'r')
#                 data = back_file[self.psr]
#             finally:
#                 back_file.close()
#                 
#             # loop over instantiations
#             for inst in data.columns:
#                 
#                 inst_number = int(n*self.background.filesize + inst)
#                 
#                 print '%i/%i ' % (inst_number, len(self.hinj)-1),
#                 
#                 
#                 # select injection psi and iota
#                 psi_inj  = random.uniform(self.pol_range[0], self.pol_range[1])
#                 iota_inj = random.uniform(self.inc_range[0], self.inc_range[1])  
#                 phi0 = random.uniform(self.phi0_range[0], self.phi0_range[1])
#                 
#                 # select instantiation
# #                 d = data[inst]
# #                 d = self.background.seed.finehet
# #                 print d
# 
#                 nlevel = 1.5e-23
#                 d = np.random.normal(scale=nlevel, size=len(self.t)) + 1j*np.random.normal(scale=nlevel, size=len(self.t))
#                 
#                 # compare sigmas
#                 plt.figure()
#                 self.sg.plot('+', label='Actual')
#                 plt.plot(self.t, np.std(d), color='r', label='Fake')
#                 plt.title('Standard deviation comparison: S6 data vs. fabricated noise')
#                 plt.legend(numpoints=1)
#                 plt.xlabel('GPS time')
#                 plt.ylabel('$\sigma$')
#                 plt.savefig('files/plots/sigmacomp.png', bbox_inches='tight')
#                 plt.close()
#                 
#                 # search
#                 h = self.hinj[inst_number]
#                 
#                 for d_inj in self.dinj:
#                     print 'I! d = %(d_inj)f' % locals()
#                     print 'gaussian noise!'
#                 
#                     # inject signal
#                     d += h * self.injection.simulate(d_inj, psi_inj, iota_inj, phase=0.0)
#                 
#                     # search
#                     res = []
#                     for delta in self.dsrch:
#                         print delta,
#                         s = signal.simulate(delta, psi, iota, phase=0.0)  
# 
#                         res += [abs(np.vdot(s, d))]
#                         
#                     self.results[d_inj][h] = self.dsrch(np.argmax(res))
#                         
#                 # save
#                 try:
#                    self.results.to_pickle('files/analysis/results/distribution')
#                 except:
#                     print "Couldn't save!"


def chi(A, b):
    # chi2 minimization through svd decomposition
    svd = np.linalg.svd(A, full_matrices=False)

    U = pd.DataFrame(svd[0], columns=A.columns, index=A.index)
    W = pd.DataFrame(np.diag(1./svd[1]), index=A.columns, columns=A.columns)
    V = pd.DataFrame(svd[2], index=A.columns, columns=A.columns)

    cov = V.T.dot(W**2).dot(V)  # covariance matrix

    VtW = V.T.dot(W)
    # need to make U complex before dotting with b
    Utb = (U + 0j).mul(b, axis=0).sum(axis=0)
    a = VtW.dot(Utb.T)          # results
    
    return a, cov

    

## SEARCH
class Sigma(object):
    def __init__(self, detector, psr, data, justload=False):
        self.detector = detector
        self.psr = psr
        
        self.data =  data
        
        self.dir = paths.sigma + '/' + self.detector + '/'
        self.name = 'segsigma_' + self.psr + '_' + self.detector
        self.path = self.dir + self.name
        
        self.justload = justload
        
        self.get()
        
    def create(self):
        '''
        Splits data into day-long segments and returns their standard deviation.
        '''
        
        data  = self.data
        
        # Check orientation    
        t = data.index
        interval_length= sd.ss
        print 'Taking std over %f second-long intervals:' % interval_length,

        # Slice up data into day-long bins and get groupby stats (see Ch 9 of Python for Data Analysis).
        bins = np.arange(t[0]-interval_length, t[-1]+interval_length, interval_length)
        slices = pd.cut(t, bins, right=False)
        print 'sliced,',

        def getsigma(group):
    #         s = np.std(group)
            g = np.array(group.tolist())
            s = np.std(g)
            return s
            #return group.std(ddof=0) # this is pd unbiased 1/(n-1), should use np.std 1/n?
        
        print 'std taken,',
        grouped = data.groupby(slices) # groups by bin
        sigmagroups= grouped.apply(getsigma) # gets std for each bin
        print 'grouped,',

        # Create standard deviation time series 
        s = [sigmagroups.ix[slices.labels[t_index]] for t_index in range(0,len(t)) ]
        self.std = pd.Series(s, index=t)
        print 'done.'
    
    
    def get(self):
        print 'Retrieving segment standard deviation...' % locals(),
        try:
            s = pd.HDFStore(self.path)
            try:
                self.std = s[self.psr]
                
                # check times coincide
                if not self.justload:
                    if not set(self.std.index)==set(self.data.index):
                        self.create()
                        # save
                        s.close()
                        s = pd.HDFStore(self.path, 'w')
                        s[self.psr] = self.std
                    
            except KeyError:
                print 'PSR not in file.',
                self.create()
                # save
                s[self.psr] = self.std
                
        except IOError:
            print 'Creating std directory.',
            os.makedirs(self.dir)
            self.create()
            # save
            s = pd.HDFStore(self.path, 'w')
            s[self.psr] = self.std
            
        finally:
            s.close()
        
        print 'Sigma is ready.'
        
    def plot(self, extra_name=''):
        
        self.std.plot(style='+')
        plt.title('Daily standard deviation for ' + self.detector + ' ' + self.psr + ' data ' + extra_name)
        plt.xlabel('GPS time (s)')
        plt.ylabel('$\sigma$')
        
        # save
        save_dir = paths.plots + '/' + self.detector + '/sigma/'
        save_name = self.name + extra_name + '.png'
        try:
            plt.savefig(save_dir + save_name, bbox_inches='tight')
        except IOError:
            os.makedirs(save_dir)
            plt.savefig(save_dir + save_name, bbox_inches='tight')