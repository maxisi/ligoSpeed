import pandas as pd
import numpy as np
import random
import csv
import os
import sys

import paths
import sidereal as sd
import templates

reload(templates)

class EphemerisD(object):

    def __init__(self, id, time, prename='earth00-19-DE', extra_name=''):
        self.id = str(id)                                   # kind of DE
        self.originalName = prename + self.id + '.dat'      # name of DE file
        self.step = 14400.0                                 # time step in originalDE
        self.name = 'DE' + self.id                          # name of imported DE
        self.extra_name = extra_name  # e.g. PSRid          # adds name to interpolation
        self.t = np.array([int(x) for x in time])
        self.nentries = len(time)
    
        self.load()                                         # load geocenter locations
        
    def get(self):
        '''
        Load Earth ephemerides DE(id) data file into DataFrame.
        '''
    
        lines_per_entry = 4
        lines_in_header = 22
    
        data_and_header = list(csv.reader(open(paths.eph + self.originalName)))
        lines = data_and_header[lines_in_header + 1:]
    
        nlines = len(lines)
        nentries = nlines/lines_per_entry
    
        x = {}; y = {}; z = {};
        vx = {}; vy = {}; vz = {};
        ax = {}; ay = {}; az = {};
    
        for i in range(nentries):
            row = {n: lines[lines_per_entry*i + n][0].split('\t') for n in range(lines_per_entry)}
            # row[0]=['', '629856000.0000000', '-0.865732865187040', '449.070722945508976']
        
            # time
            t = int(float(row[0][1]))            # 629856000.0
        
            # radius
            x[t] = float(row[0][2])              # -0.86573286518704
            y[t] = float(row[0][3])              # 449.070722945508976
            z[t] = float(row[1][1])              # 194.808873897521295
        
            # velocity
            vx[t] = float(row[1][2])
            vy[t] = float(row[1][3])
            vz[t] = float(row[2][1])

            # acceleration
            ax[t] = float(row[2][2])
            ay[t] = float(row[2][3])
            az[t] = float(row[2][1])

        # frame
        r = pd.DataFrame([x, y, z], index = ['x', 'y', 'z'])
        v = pd.DataFrame([vx, vy, vz], index = ['x', 'y', 'z'])
        a = pd.DataFrame([ax, ay, az], index = ['x', 'y', 'z'])
    
        # store
        try:
            f = pd.HDFStore(paths.eph_local + self.originalName, 'w')
            f['r'] = r
            f['v'] = v
            f['a'] = a
            print 'Imported DE %s.' % self.id
        finally:
            f.close()


    def produce(self):
        '''
        Loads Earth ephemeris data and interpolates to obtain Earth SSB location for time.
        '''
    
        # Try to load DE in HDF5. Generate file if necessary
        try:
            f = pd.HDFStore(paths.eph_local + self.originalName, 'r')
        except IOError:
            self.get()
            f = pd.HDFStore(paths.eph_local + self.originalName, 'r')
        finally:
            rDE = f['r']
            vDE = f['v']
            aDE = f['a']
            f.close()

        # Ephemeris time
        tDE = np.array(rDE.columns)
        
        # Indices of closest DE times to arrival times
        tDEi = np.floor( (self.t - tDE[0])/self.step )
        
        # Time selection from DE
        tFilter = tDE[tDEi.astype(int)]
        tF = {self.t[i]: tFilter[i] for i in range(self.nentries)}
        
        # Get difference between the two
        dt = self.t - tFilter
        dT = {self.t[i]: dt[i] for i in range(self.nentries)}

        dt2 = (dt**2)/2.
        dT2 = {self.t[i]: dt2[i] for i in range(self.nentries)}
        
        # Interpolate
        self.r = pd.DataFrame({ti: rDE[tF[ti]] + vDE[tF[ti]]*dT[ti] + aDE[tF[ti]]*dT2[ti] for ti in self.t})
        
        try:
            f = pd.HDFStore(paths.eph_local + self.name + self.extra_name, 'w')
            f['r'] = self.r
            print 'Interpolation successful. ',
        finally:
            f.close()
            
            
    def fileload(self):
        '''
        Tries to load interpolated files. Calls interpolation if this fails
        '''
        try:
            f = pd.HDFStore(paths.eph_local + self.name + self.extra_name, 'r')
            self.r = f['r']
            f.close()
        except IOError:
            # Interpolate positions (not loaded and failed to get file)
            print 'Interpolating ephemeris. ',
            self.produce()
            
            
    def load(self):
        '''
        Loads interpolated SSB geocenter locations.
        '''
        # Check if the positions are loaded already
        try:
            self.r
        except AttributeError:
            # Check if positions are stored in file
            self.fileload()
            
        # Check data type
        try:
            tDE = self.r.columns
        except AttributeError:
            self.fileload()
            tDE = self.r.columns

        try:  
            if set(tDE)==set(self.t):
                print 'All times present. ',
            else:
                print 'Interpolating ephemeris. ',
                self.produce()
        except TypeError:
                print 'Interpolating ephemeris (TypeError). ',
                self.produce()

        print 'Geocenter positions ready.'


class EphemerisT(object):

    def __init__(self, time, prename='te405_2000-2019.dat', extra_name=''):
        self.originalName = prename                         # name of DE file
        self.name = 'TE405'                                 # name of imported DE
        self.step = 14400.0                                 # time step in originalDE
        self.extra_name = extra_name  # e.g. PSRid          # adds name to interpolation
        self.t = np.array([int(x) for x in time])
        self.nentries = len(time)
    
        self.load()                                         # load geocenter locations
        
    def get(self):
        '''
        Load time ephemerides TE405 data file into DataFrame.
        '''
        
        teF = pd.read_table(paths.eph + originalName, sep='\n', skiprows=17, header=None, names='T')
        te = teF['T']
        te.index = np.linspace(630720000, 630720000+14400*43890, 43890)

        # store
        try:
            f = pd.HDFStore(paths.eph_local + self.originalName, 'w')
            f['t'] = te
            print 'Imported TE405.'
        finally:
            f.close()


    def produce(self):
        '''
        Loads Earth ephemeris data and interpolates to obtain Earth SSB location for time.
        '''
            
        # Try to load TE from HDF5. Generate file if necessary
        try:
            print '1'
            f = pd.HDFStore(paths.eph_local + self.originalName, 'r')
        except IOError:
            print '2'
            self.get()
            f = pd.HDFStore(paths.eph_local + self.originalName, 'r')
        finally:
            print '3'
            tE = f['t']
            f.close()
        
        # Indices of closest TE times to arrival times
        tEi = np.floor( (self.t - tE[0])/self.step )
        
        # Time selection from TE
        self.tE = tE[tEi.astype(int)]
        
        try:
            f = pd.HDFStore(paths.eph_local + self.name + self.extra_name, 'w')
            f['t'] = self.tE
            print 'Interpolation successful. ',
        finally:
            f.close()
            
            
    def fileload(self):
        '''
        Tries to load interpolated files. Calls interpolation if this fails
        '''
        try:
            f = pd.HDFStore(paths.eph_local + self.originalName + self.extra_name, 'r')
            self.tE = f['t']
            f.close()
        except IOError:
            # Interpolate positions (not loaded and failed to get file)
            print 'Interpolating ephemeris. ',
            self.produce()
            
            
    def load(self):
        '''
        Loads interpolated SSB geocenter locations.
        '''
        # Check if the positions are loaded already
        try:
            self.tE
        except AttributeError:
            # Check if positions are stored in file
            self.fileload()
            
        # Check data type
        try:
            tDE = self.r.columns
        except AttributeError:
            self.fileload()
            tDE = self.r.columns
            
        if all(tDE==self.t):
            print 'All times present. ',
        else:
            print 'Interpolating ephemeris. ',
            self.produce()
        
        print 'Geocenter positions ready.'


class System(object):
    def __init__(self, detector, psr):
        self.detector = detector
        self.det = sd.detnames(detector)
        self.psr = psr
        
        self.src = templates.Source(psr)
        self.obs = templates.Detector(detector)
        
    def getroemer(self, t, c=sd.c):
        
        coords = ['x', 'y', 'z']
        
        # vector from geocenter to detector
        self.obs.t = np.array(t).astype(int)
        self.obs.loadVectors()
        r_d = self.obs.dz
        
        # unit vector pointing to source (opposite of source-SSB):
        # building following Edwards et al. TEMPO2 (14)
        n0 = - self.src.wz
        
        aL = [
                -np.sin(self.src.param['RAS']),
                np.cos(self.src.param['DEC']),
                0.
            ]
        
        dL = [
                -np.cos(self.src.param['RAS'])*np.sin(self.src.param['DEC']),
                -np.sin(self.src.param['RAS'])*np.sin(self.src.param['DEC']),
                np.cos(self.src.param['DEC'])
            ]

        a = pd.Series(aL, index=coords)
        d = pd.Series(dL, index=coords)

       
        mu_perp = self.src.param['PMRAS']*a + self.src.param['PMDEC']*d
        
        dt = pd.Series(self.obs.t - self.src.param['POSEPOCH'], index=self.obs.t)
        
        speed_term = pd.DataFrame(np.outer(mu_perp, dt), index=coords, columns=self.obs.t)
        
        accel_term_A = - abs(mu_perp.dot(mu_perp))**2 * n0 /2. 
        accel_term = pd.DataFrame(np.outer(accel_term_A, dt**2), index=coords, columns=self.obs.t)
        
        n = (speed_term + accel_term).add(n0, axis='index').T
        
        
        # vector from SSB to geocenter
        if self.psr == 'J0534+2200':
            ephID = 200
        else:
            ephID = 405
        
        eph = EphemerisD(ephID, self.obs.t)
        
        # vector from SSB to detector
        r = r_d + eph.r.T
        
        # dot product
        rn = (r*n).sum(axis=1)
        
        self.roemer = rn/c

    def einstein(self, tgps, tdb=False):
        '''
        Returns TCB - TT at a given arrival time. If tdb=True, returns TDB - TT instead
        i.e. a scaled version of the TCB frame in which the mean drift relative to TT is
        divided out (TEMPO2 doc).

        Stolen from XLALBarycenter
        '''
    
        if tdb:  
            jedtdt = -7300.5e0 + (tgps + 51.184)/8.64e4; 
            jt=jedtdt/3.6525e5; # converting to TEMPO expansion param Julian millenium, NOT Julian century

            dE = 1.e-6*(
            1656.674564e0 * np.sin(6283.075849991e0*jt + 6.240054195e0 )+ 
            22.417471e0 * np.sin(5753.384884897e0*jt + 4.296977442e0 ) + 
            13.839792e0 * np.sin(12566.151699983e0*jt + 6.196904410e0 ) + 
            4.770086e0 * np.sin(529.690965095e0*jt + 0.444401603e0 ) + 
            4.676740e0 * np.sin(6069.776754553e0 *jt + 4.021195093e0 ) + 
            2.256707e0 * np.sin(213.299095438e0 *jt + 5.543113262e0 ) + 
            1.694205e0 * np.sin(-3.523118349e0 *jt + 5.025132748e0 ) + 
            1.554905e0 * np.sin(77713.771467920e0 *jt + 5.198467090e0 ) + 
            1.276839e0 * np.sin(7860.419392439e0 *jt + 5.988822341e0 ) + 
            1.193379e0 * np.sin(5223.693919802e0 *jt + 3.649823730e0 ) + 
            1.115322e0 * np.sin(3930.209696220e0 *jt + 1.422745069e0 ) + 
            0.794185e0 * np.sin(11506.769769794e0 *jt + 2.322313077e0 ) + 
            0.447061e0 * np.sin(26.298319800e0 *jt + 3.615796498e0 ) + 
            0.435206e0 * np.sin(-398.149003408e0 *jt + 4.349338347e0 ) + 
            0.600309e0 * np.sin(1577.343542448e0 *jt + 2.678271909e0 ) + 
            0.496817e0 * np.sin(6208.294251424e0 *jt + 5.696701824e0 ) + 
            0.486306e0 * np.sin(5884.926846583e0 *jt + 0.520007179e0 ) + 
            0.432392e0 * np.sin(74.781598567e0 *jt + 2.435898309e0 ) + 
            0.468597e0 * np.sin(6244.942814354e0 *jt + 5.866398759e0 ) + 
            0.375510e0 * np.sin(5507.553238667e0 *jt + 4.103476804e0 ) 
            )

            dE += 1.e-6*(
            0.243085 * np.sin(-775.522611324 *jt + 3.651837925 )   +
            0.173435 * np.sin(18849.227549974 *jt + 6.153743485 )   +
            0.230685 * np.sin(5856.477659115 *jt + 4.773852582 )   +
            0.203747 * np.sin(12036.460734888 *jt + 4.333987818 )   +
            0.143935 * np.sin(-796.298006816 *jt + 5.957517795 )   +
            0.159080 * np.sin(10977.078804699 *jt + 1.890075226 )   +
            0.119979 * np.sin(38.133035638 *jt + 4.551585768 )   +
            0.118971 * np.sin(5486.777843175 *jt + 1.914547226 )   +
            0.116120 * np.sin(1059.381930189 *jt + 0.873504123 )   +
            0.137927 * np.sin(11790.629088659 *jt + 1.135934669 )   +
            0.098358 * np.sin(2544.314419883 *jt + 0.092793886 )   +
            0.101868 * np.sin(-5573.142801634 *jt + 5.984503847 )   +
            0.080164 * np.sin(206.185548437 *jt + 2.095377709 )   +
            0.079645 * np.sin(4694.002954708 *jt + 2.949233637 )   +
            0.062617 * np.sin(20.775395492 *jt + 2.654394814 )   +
            0.075019 * np.sin(2942.463423292 *jt + 4.980931759 )   +
            0.064397 * np.sin(5746.271337896 *jt + 1.280308748 )   +
            0.063814 * np.sin(5760.498431898 *jt + 4.167901731 )   +
            0.048042 * np.sin(2146.165416475 *jt + 1.495846011 )   +
            0.048373 * np.sin(155.420399434 *jt + 2.251573730 )
            )
        
        else:
            eph = EphemerisT(tgps)
            tE = eph.Te
            # assuming tgps is in fact the time ephemeris
            deltaT = np.array(list(tE))
            tgps = np.array(tgps.index)
    
            IFTE_KM1 = 1.55051979176e-8
            IFTE_K = 1. + IFTE_KM1
            IFTE_MJD0 = 43144.0003725
            IFTE_TEPH0 = -65.564518e-6
            IFTE_LC = 1.48082686742e-8
        
            scorr = IFTE_K # assuming TCB input

            mjdtt = 44244. + (tgps + 51.184)/86400.
    
        #     deltaT = 0.0001332282928461 #random value taken from real eph. just to get magnitude.
    
            correctionTT_Teph = IFTE_TEPH0 + deltaT / (1.0-IFTE_LC);
    
            dE = IFTE_KM1 * (mjdtt - IFTE_MJD0)*86400.0 + IFTE_K * (correctionTT_Teph - IFTE_TEPH0);
    #     
        return dE
     
     
    
def fakedata(ndays, scale=10**-21, t0=630720013):
    t = np.arange(t0, t0 + ndays*sd.ss, 1000*sd.periodLIGO).astype(int)
    d = [random.random()*scale for x in t]
    data = pd.Series(d, index=t)
    return data
    

class DataReduction(object):
    def __init__(self, detector, psr, ndays, c_gw=sd.c):
        # Detector and source
        self.detector = detector
        self.psr = psr
        self.syst = System(detector, psr)
        
        # noise
        print 'Generating noise.'
        self.noise = fakedata(ndays)
        
        # delays
        self.t = self.noise.index
        print 'Computing Roemer delay.'
        self.syst.getroemer(self.t, c=c_gw)
        
        # signal
        print 'Composing signal.'
        s = templates.Signal(self.detector, self.psr, 'GR', 'p', self.t)
        self.h = s.simulate(self.syst.src.param['POL'], self.syst.src.param['INC'])
        
        print 'Adding signal to noise.'
        self.s = self.h + self.noise