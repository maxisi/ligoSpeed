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
                self.imp()
            finally:
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
    def __init__(self, detector, psr, dinj=[], hinj=[], pdif_s=None, kind='GR', pdif=None):
        # system
        self.detector = detector
        self.psr = psr
        
        # search
        self.methods = ['GR']
        
        # injection
        self.dinj = dinj
        self.hinj = hinj
        
        self.kind = kind
        
        self.pdif = pdif
        self.pdif_s = pdif_s
        
        # parameters
        self.parameters = ['h', 's', 'drec', 'sphi']
        
        # containers
        for p in self.parameters:
            setattr(self, p, pd.DataFrame(columns = dinj, index = hinj))
                
        # saving
        self.dir = paths.results + self.detector + '/' + self.psr + '/' 
        self.name = self.psr + '_' + self.detector + '_' + self.kind + '_' + sd.phase2(pdif)
        self.path = self.dir + self.name
        
        self.issaved =  False
        
    def save(self):
        
        self.h.index = self.hinj
        self.s.index = self.hinj
               
        try:
            os.makedirs(self.dir)
        except OSError:
            pass
            
        try:
            f = pd.HDFStore(self.path, 'w')
            
            for p in self.parameters:
                f[p] = getattr(self, p)

            f['stats']= self.stats
            
        finally:
            f.close()
            
        self.issaved = True
        
        
    def load(self):
        try:
            f = pd.HDFStore(self.path, 'r')
            for p in self.parameters:
                setattr(self, p, f[p])
            self.s = f['s']
        finally:
            f.close() 


    def plots(self, pltType, extra_name=''):
        
        header = self.kind + sd.phase2(self.pdif) + ' injections on ' + self.detector + ' data for ' + self.psr + ' ' + extra_name
          
        getattr(psrplot, pltType)(hinj=self.h.index, hrec=self.h, s=self.s, methods=self.methods)
        
        plt.title(header)
        
        pltdir = paths.plots + self.detector + '/' + self.kind + '/' + pltType + '/'
        pltname = self.detector + '_' + self.kind + '_' + sd.phase2(self.pdif) + '_' + pltType + extra_name
        save_to = pltdir + pltname
        
        try:
            os.makedirs(pltdir)
        except OSError:
            pass
            
        plt.savefig(save_to, bbox_inches='tight')
        plt.close()
        
        print 'Plot saved to:\n %(save_to)s' % locals()
    
        
    def getstats(self, plot=False, store=True):
        
        lins = self.s.applymap(math.sqrt)
        
        for m in self.methods:
            self.stats[m]['min inj det'] = psrplot.min_det_h(lins[m])
            self.stats[m]['lin s slope'] = psrplot.lin_fit(lins[m])(1)
            self.stats[m]['lin s noise'] = psrplot.noise_line(lins[m])(1)
            self.stats[m]['lin s inter'] = psrplot.fit_intersect_noise(lins[m])
            self.stats[m]['h rec noise'] = psrplot.noise_line(self.h[m])(1)
            self.stats[m]['h rec slope'] = psrplot.lin_fit(self.h[m])(1)
            self.stats[m]['h rec inter'] = psrplot.fit_intersect_noise(self.h[m])


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



class InjSearch(object):
    
    def __init__(self, detector, psr, nfreq, pdif, nd, sweep=101, rangeparam=[], dinjrange=[1., 1.01], dsrchrange=[1., 1.01], hinjrange=[1.0E-27, 1.0E-23], filesize=100):
        # system info
        self.detector = detector
        self.psr = psr
                
        # data info
        frange = [1.0e-7, 1.0e-5]
        self.freq = np.linspace(frange[0], frange[1], nfreq)
        print 'Getting background.'
        self.background = Background(detector, psr, self.freq, filesize)
        self.background.get()
        
        self.t = self.background.seed.finehet.index
        
        sigma = Sigma(self.detector, self.psr, self.background.seed.finehet)
        self.sg = sigma.std
        
        # injection strengths
        self.hinj = np.linspace(hinjrange[0], hinjrange[1], nfreq)
        
        # injection deltas
        self.dinj = np.linspace(dinjrange[0], dinjrange[1], nd)
        
        # injection kind   
        self.pdif = pdif
        self.injkind = 'GR'
        
        # set up injection
        self.injection = templates.Signal(detector, psr, pdif, self.t)
        self.setranges(rangeparam)
        
        # search
        self.dsrch = np.linspace(dsrchrange[0], dsrchrange[1], sweep)
        
        
    def setranges(self, rangeparam):
        src = self.injection.response.src
        
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
    
    
    def analyze(self):

        print 'Analyzing %d files.' % self.background.nsets
    
        # search info
        search = templates.Signal(self.detector, self.psr, self.pdif, self.t)

        # results
        self.results = Results(self.detector, self.psr, dinj=self.dinj, hinj=self.hinj, pdif=self.pdif)
            
        # loop over files
        for n in range(self.background.nsets):
            
            try:
                back_file = pd.HDFStore(self.background.path + str(n), 'r')
                data = back_file[self.psr]
            finally:
                back_file.close()
                
            # loop over instantiations
            for inst in data.columns:
                
                inst_number = int(n*self.background.filesize + inst)
                
                print '%i/%i ' % (inst_number, len(self.hinj)-1),
                
                # select search psi, iota and phi0
                psi  = random.uniform(self.pol_range[0], self.pol_range[1])
                iota = random.uniform(self.inc_range[0], self.inc_range[1])
                
                # select injection psi and iota
                psi_inj  = random.uniform(self.pol_range[0], self.pol_range[1])
                iota_inj = random.uniform(self.inc_range[0], self.inc_range[1])  
                phi0 = random.uniform(self.phi0_range[0], self.phi0_range[1])
                
                # phi0 range to sweep over        
#                 phi0_srch = np.linspace(self.phi0_range[0], self.phi0_range[1], 5)
                
                # select instantiation
                d = data[inst]

                h = self.hinj[inst_number]

                if h != 0:
                
                    for d_inj in self.dinj:
                    
                        print 'I! d = %(d_inj)f\t%(psi_inj)f %(iota_inj)f %(phi0)f' % locals()
                        # print 'no noise!'
                        
                        # inject signal
                        d += h * self.injection.simulate(d_inj, psi_inj, iota_inj, phase=phi0)
                        
                        # search
                        res = pd.Series(index=self.dsrch)
                        for delta in self.dsrch:
                        
                            s = search.simulate(delta, psi_inj, iota_inj, phase=phi0)  
                            
                            res[delta] = abs(np.vdot(d, s))
                            
                        self.results.drec[d_inj][h] = res.index[np.argmax(res)]

        ## Save
        self.results.save()


class Distribution(InjSearch):
    
    def __init__(self):
        super(Distribution, self).__init__('H1', 'J0534+2200', 10, 0., 50)
        
    def analysis(self):

        print 'Getting matched filter distribution.'
    
        # search info
        search = templates.Signal(self.detector, self.psr, self.pdif, self.t)

        # results
        self.results = pf.DataFrame(columns=self.dinj, index=self.dsrch)
            
        # loop over files
        for n in range(self.background.nsets):
            
            try:
                back_file = pd.HDFStore(self.background.path + str(n), 'r')
                data = back_file[self.psr]
            finally:
                back_file.close()
                
            # loop over instantiations
            for inst in data.columns:
                
                inst_number = int(n*self.background.filesize + inst)
                
                print '%i/%i ' % (inst_number, len(self.hinj)-1),
                
                
                # select injection psi and iota
                psi_inj  = random.uniform(self.pol_range[0], self.pol_range[1])
                iota_inj = random.uniform(self.inc_range[0], self.inc_range[1])  
                phi0 = random.uniform(self.phi0_range[0], self.phi0_range[1])
                
                # phi0 range to sweep over        
#                 phi0_srch = np.linspace(self.phi0_range[0], self.phi0_range[1], 5)
                
                # select instantiation
                d = data[inst]

                h = self.hinj[inst_number]

                if h == max(self.hinj):
                
                    for d_inj in self.dinj:
                    
                        print 'I! d = %(d_inj)f\t%(psi_inj)f %(iota_inj)f %(phi0)f' % locals()
                        # print 'no noise!'
                        
                        # inject signal
                        d += h * self.injection.simulate(d_inj, psi_inj, iota_inj, phase=phi0)
                        
                        # search
                        res = pd.Series(index=self.dsrch)
                        for delta in self.dsrch:
                        
                            s = search.simulate(delta, psi_inj, iota_inj, phase=phi0)  
                            
                            self.results[d_inj][delta] = abs(np.vdot(d, s))
    

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