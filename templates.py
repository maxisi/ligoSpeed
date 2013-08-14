import math
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from urllib2 import urlopen
from collections import namedtuple
import csv

import sys, os
from sys import argv, exit
import copy
from time import time

import sidereal as sd
import paths

reload(sd)
    
#  vector container to make it easy to move vectors around
Vectors = namedtuple('Vectors', ['dx', 'dy', 'wx', 'wy', 'wz'])


## SOURCE
def currentCatalogue():
    '''
    Returns contents of current pulsar catalogue. Mainly for maintenance reasons.
    '''
    try:
        return pd.read_pickle(paths.psrcat)
    except:
        print 'No pulsar catalogue found.'
        exit()
        

class Source(object):
    '''
    Assuming single pulsar.
    '''
    
    def __init__(self, psr, loadvec=True):
            
        self.psr = psr
        self.npsrs = 1
        self.path = paths.vectors + 'srcVec_' + psr
                     
        # If necessary, rebuild catalogue; otherwise, just load catalogue.
        try:
            f = open(paths.textfromATNF, 'r')
            pd.read_pickle(paths.psrcat)
        except IOError:
            self.build_catalogue()
            f = open(paths.textfromATNF, 'r')
        finally:
            f_text = f.read()
            f.close()
        
        if self.psr not in f_text:
            self.build_catalogue()

        psrcat = pd.read_pickle(paths.psrcat)
        self.param = psrcat.ix[self.psr]
        
        if loadvec:
            self.loadVectors()
        

    def build_catalogue(self, extrapsrs=[]):
        '''
        Gets location parameters for pulsars in input list from ATNF online catalogue.
        Creates and pickles corresponding pandas DataFrame 'pulsar_parameters'.
        '''
        psrs = [self.psr] + extrapsrs
        
        def atnfurl(psr_list):
            pre = 'http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.47&JName=JName&RaJ=RaJ&DecJ=DecJ&PMRA=PMRA&PMDec=PMDec&PosEpoch=PosEpoch&F0=F0&F1=F1&F2=F2&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&'
            post = '&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+errors&no_value=*&nohead=nohead&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=64&table_bottom.y=24'
            names= 'pulsar_names='
    
            for psr in psr_list:
                names+=psr.replace('+', '%2B')
        
                if psr != psr_list[-1]:
                    names+='%0D%0A'
            url = "%(pre)s%(names)s%(post)s" % locals()
    
            if len(url)>2000:
                print 'WARNING! URL %d characters!' % len(url)
    
            return url

        # Get data
        url = atnfurl(psrs)
        soup = BeautifulSoup(urlopen(url)) # get webpage and parse
        text = str(soup.pre.string)

        # Write txt file 
        f = open(paths.textfromATNF, 'w')
        f.write(text)
        f.close()

        # Load and create DataFrame
        psr = pd.read_table(paths.textfromATNF, sep='\s+', comment='*', names=sd.paramNames, header=None, skiprows=1, index_col=1)
        psrcat=psr.drop("#", axis=1)
    
        # Format
        formatRAS  = lambda x: sd.hms_rad(x)
        formatRASe = lambda x: sd.hms_rad(0., 0., x)
        formatDEC  = lambda y: np.radians(sd.dms_deg(y))
        formatDECe = lambda x: np.radians(sd.dms_deg(0., 0., x))
        
        formatPM    = lambda x: sd.masyr_rads(x)
        formatEPOCH = lambda x: sd.mjd_gps(x)
        
        psrcat['RAS'] = psrcat['RAS'].map(formatRAS)                    # hms -> rad
        psrcat['RAS error'] = psrcat['RAS error'].map(formatRASe)
        psrcat['DEC'] = psrcat['DEC'].map(formatDEC)                    # dms -> rad
        psrcat['DEC error'] = psrcat['DEC error'].map(formatDECe) 
        
        psrcat['PMRAS'] = psrcat['PMRAS'].map(formatPM)                 # mas/yr -> rad/s
        psrcat['PMRAS error'] = psrcat['PMRAS error'].map(formatPM)
        psrcat['PMDEC'] = psrcat['PMDEC'].map(formatPM)
        psrcat['PMDEC error'] = psrcat['PMDEC error'].map(formatPM)
        
        psrcat['POSEPOCH'] = psrcat['POSEPOCH'].map(formatEPOCH)        # MJD -> GPS

    
        #Check extra parameters
        extra = pd.read_table(paths.psrextra, sep=',', header=None, index_col=[0], names=sd.extraParamNames)
        psrcatComplete = pd.merge(psrcat,extra, left_index=True, right_index=True, how='outer')
    
        psrcatComplete['POL error']=psrcatComplete['POL error'].fillna(value=np.pi/4.)
        psrcatComplete['INC error']=psrcatComplete['INC error'].fillna(value=np.pi/4.)

    
        psrcatComplete.fillna(value=0, inplace=True)
        psrcatComplete.save(paths.psrcat)
     
     
    def loadVectors(self):
        '''
        Loads detector arm vectors from file.
        '''
        try:
            file = pd.HDFStore(self.path, 'r')
            self.wx = file['wx']
            self.wy = file['wy']
            self.wz = file['wz']
            file.close()
        except IOError:
            self.createVectors()
           
            
    def createVectors(self):
        # Return wave vectors for all sources listed.
        
        north = np.array([0, 0, 1])
        
        # take source location vector components in celestial coordinates and invert direction multiplying by -1 to get wave vector wz
        wz = [-math.cos(self.param['DEC'])*math.cos(self.param['RAS']), -math.cos(self.param['DEC'])*math.sin(self.param['RAS']), -math.sin(self.param['DEC'])] 
        self.wz = pd.Series(wz, name=self.psr, index=['x', 'y', 'z'])
        
        
        wy = np.cross(wz, north) 
        wy /= np.sqrt(np.sum(wy ** 2))
        self.wy = pd.Series(wy, name=self.psr, index=['x','y','z'])

        wx = np.cross(wy, wz)
        wx /= np.sqrt(np.sum(wx ** 2))
        self.wx = pd.Series(wx, name=self.psr, index=['x','y','z'])
                
        try:
            f = pd.HDFStore(self.path, 'w')
            f['wx'] = self.wx
            f['wy'] = self.wy
            f['wz'] = self.wz
        finally:
            f.close()
         
###         

# DETECTOR
class Detector(object):
    
    def __init__(self, d, t=[]):
        self.id = d  
        self.name = sd.detnames(d)
        self.param = sd.detectors[self.name]
        self.t = np.array(t).astype(int)
        self.nentries = len(t)
        self.path = paths.vectors + 'detVec_' + self.name
        
        if self.t!=[]:
            self.loadVectors()
        

    def fileload(self):
        '''
        Loads detector arm vectors from file.
        '''
        try:
            file = pd.HDFStore(self.path, 'r')
            self.dx = file['dx']
            self.dy = file['dy']
            self.dz = file['dz']
            file.close() 
        except IOError:
            self.createVectors()


    def loadVectors(self):
        '''
        Loads detector arm vectors if necessary.
        '''
        # Check if vectors are stored in file'
        self.fileload()
            
        # Check data type'
        try:
            self.dx.columns
        except AttributeError:
            self.fileload()
        
#         try:
        if set(self.dx.index)==set(self.t):
            # All times present.
            pass
        else:
            self.createVectors()
#         except TypeError:
#             self.createVectors()

    def createVectors(self):
        '''
        Returns arm vectors in Cartesian sidereal coordinates.
        '''
        print 'Creating detector vectors.'
        northPole = pd.Series(np.array([0, 0, 1]), index=['x', 'y', 'z'])    # Earth center to North pole
        
        lat = self.param.ix['lat']
        lon = self.param.ix['lon']
        x_east = self.param.ix['x_east']
        arm_ang = self.param.ix['arm_ang']
        
        coords = ['x', 'y', 'z']
        
        t = np.array(self.t).astype(int)

        length = self.nentries

        # Angle between detector and Aries (vernal equinox) at time t
        # fiducial GPS time t0=630763213 (12hUT1 1/1/2000, JD245154).
        # See http://aa.usno.navy.mil/faq/docs/GAST.php
        offset = 67310.5484088 * sd.w   # Aries-Greenwich angle at fiducial time (GMST)
        lmst = offset + sd.w*(t-630763213) + lon # (LMST) rows: t, columns: det
        th = pd.Series(lmst, index=t)
        
        start = time()
        z = {
            'x' : np.cos(lat)*np.cos(th),
            'y' : np.cos(lat)*np.sin(th),
            'z' : pd.Series([math.sin(lat)]*len(t), index=t)
            }  # [[x0, ...], [y0, ...], [z0, ...]]
          
        zenith = pd.DataFrame(z) # norm 1 already
        
        z2 = [zenith['x'], zenith['y'], zenith['z']]

        localEast = np.cross(northPole, z2 , axisb=0)
        localNorth = np.cross(z2, localEast, axisa=0)   
        
        xArm = np.cos(x_east)*localEast + np.sin(x_east)*localNorth
        xArm /= np.sqrt(np.sum(xArm ** 2., axis=1))[..., None]

        perp_xz = np.cross(zenith, xArm)
        yArm = xArm*np.cos(arm_ang) + perp_xz*np.sin(arm_ang) # equals perp_xz when angle between arms is 90deg

        self.dx = pd.DataFrame(xArm, index=t, columns= ['x', 'y', 'z'])
        self.dy = pd.DataFrame(yArm, index=t, columns= ['x', 'y', 'z'])
        self.dz = sd.rE * zenith

        try:
            f = pd.HDFStore(self.path, 'w')
            f['dx'] = self.dx
            f['dy'] = self.dy
            f['dz'] = self.dz
        finally:
            f.close()


###

# ANTENNA PATTERNS
class Polarizations(object):
    '''
    Produces detector response for different polarization given input detector and
    source vectors. Vectors must be in namedtuple form as defined above.
    Effectively computes dyadic product between wave and detector tensors.
    '''
    
    def __init__(self, vectors):
        # assuming vectors is of the Vectors container kind
        self.vec = vectors
    
                    
    def product(self, polKey):
        # check all necessary vector products exist. Otherwise, create them.
        for pair in sd.polComponents[polKey]:
        
            pairName = pair[0] + pair[1]    # 'wxdx' = ('wx','dx')[0] + ('wx','dx')[1]
            
            # check if product is already defined.
            if pairName not in dir(self):
            
                # get vectors
                v0 = getattr(self.vec, pair[0])     # v0 = self.vec.wx
                v1 = getattr(self.vec, pair[1])     # v1 = self.vec.dx
                
                # dot product (mind the order! v1 MUST be detector vector to broadcast)
                setattr(self, pairName, v1.mul(v0, axis='columns').sum(axis=1))
    
    # tensor
    def plus(self):
        self.product('pl')
        pl = (self.wxdx**2 - self.wxdy**2 - self.wydx**2 + self.wydy**2)/2.
        return pl
        
    def cross(self):
        self.product('cr')
        cr = self.wxdx*self.wydx - self.wxdy*self.wydy
        return cr
    
    # vector
    def vector_x(self):
        self.product('xz')
        xz = self.wxdx*self.wzdx - self.wxdy*self.wzdy
        return xz
    
    def vector_y(self):
        self.product('yz')
        yz = self.wydx*self.wzdx - self.wydy*self.wzdy
        return yz
    
    # scalar
    def breathing(self):
        self.product('br')
        br = np.sqrt(2)*(self.wxdx**2 - self.wxdy**2 + self.wydx**2 - self.wydy**2)/2.
        return br
        # Added factor of sqrt(2) to distribute power equally among polarizations.
        # Same for longitudinal.
        
    def longitudinal(self):
        self.product('lo')
        lo = (np.sqrt(2)*(self.wzdx**2 - self.wzdy**2))/2.
        return lo
        # Modified:1/2 (based on derivation of dyadic products using tensors shown in 
        # "Gravitational wave polarizations" by Bryant Garcia.
        # The factor of 2 shouldn't be there)


class Response(object):
    '''
    Contains response of 'det' to signals from source 'psr' of kind 'kinds', over 't'.
    
    For default polarization angle, use 'get' method. This will try to recover patterns
    from file and will generate them if not found. It rotates source vectors by
    polarization angles established on file. After 'get' is called, the patterns are
    stored in the class, so there is no need to call again in same session.
    
    To input polarization angle, call 'create' method directly with argument 'psi='.
    This will ignore existing APs (if any) and will fall back on detector vectors to
    compute response (after rotation by 'psi'). Call with 'savefile=False' to prevent
    method from storing results. This routine will compute dyadic products once, so it
    can be called multiple times a session.
    
    To do: Add option to select psi randomly from range?
    '''
    
    def __init__(self, psr, det, t, kinds, loadvectors=False):
        self.det = sd.detnames(det)
        self.t = np.array(t)
        
        self.psr = psr
        self.path = paths.ap + 'ap' + psr + '_' + self.det
        
        if kinds in sd.tempNames:
            # requesting bases for preset template
            self.kinds = sd.aps[kinds]
        elif isinstance(kinds, list) and all([k in sd.tempNames for k in kinds]):
            # list of preset templates
            self.kinds = sum([sd.aps[temp] for temp in kinds], [])
        elif all([k in sd.names for k in kinds]) and not isinstance(kinds, basestring):
            # requesting bases by name
            self.kinds = kinds
        elif isinstance(kinds, basestring) and kinds in sd.names:
            self.kinds = [kinds]
        else:
            # wrong input
            print 'ERROR: %(kinds)s is not recognized as a valid basis. ap19' % locals()
            print 'Valid bases are:'
            print '\t %r' % sd.names
            print '\t %r' % sd.aps
            exit()
            
        self.hasvectors = False
        self.haspatterns = False
        
        if loadvectors:
            self.src = Source(self.psr)
            self.obs = Detector(self.det, self.t)
            self.hasvectors = True                    
                
    def create(self, savefile=True, psi=[]):
    
        # Retrieve detector vectors
        if not self.hasvectors:
            self.src = Source(self.psr)
            self.obs = Detector(self.det, self.t)
            self.hasvectors = True
            
        # Rotate source vectors
        if psi==[]:
            self.psi = self.src.param['POL']
        else:
            self.psi = psi

        wxRot = -self.src.wy*np.cos(self.psi) + self.src.wx*np.sin(self.psi)
        wyRot = self.src.wx*np.cos(self.psi) + self.src.wy*np.sin(self.psi) 
        
        # Package vectors
        vecs = Vectors(self.obs.dx, self.obs.dy, wxRot, wyRot, self.src.wz)      
           
        # Get polarizations
        pols = Polarizations(vecs)
        [setattr(self, k, getattr(pols, sd.polNames[k])() ) for k in self.kinds]
            
        # Save if requested
        if savefile:
            try:
                apF = pd.HDFStore(self.path, 'w')
                for k in self.kinds:
                    apF[k] = getattr(self, k)
            finally:
                apF.close()
                
        self.haspatterns = True
    
    def get(self):
        # Assumes no APs loaded. Otherwise, will re-write.
        
        # check disk
        try:
            apFile = pd.HDFStore(self.path, 'r')
            
        except IOError:
            # no patterns on file
            self.create()
            
        else:
            # file found
            try:
                filePols = [s.strip('/') for s in apFile.keys()]
            
                # Check file is alright:
                allPolsPresent = all([k in filePols for k in self.kinds])
            
                timeCoincides = set(self.t.astype(int)) == set(apFile[filePols[0]].index)
                # if file is empty, raises IndexError
            
                if allPolsPresent and timeCoincides:
                    [setattr(self, p, apFile[p]) for p in self.kinds]
                else:
                    apFile.close()
                    self.create()
                    
            except IndexError:
                apFile.close()
                self.create()
            
            finally:
                apFile.close()
            self.haspatterns = True        
        
    def exportAPmatlab(psr, detname, t):
        p = ap.getAP('LHO', t)
    
    #     psrdict = {
    #                 'pl':pl.T[psr].tolist(),
    #                 'cr':cr.T[psr].tolist(),
    #                 'br':br.T[psr].tolist(),
    #                 'lo':lo.T[psr].tolist(),
    #                 'xz':xz.T[psr].tolist(),
    #                 'yz':yz.T[psr].tolist()
    #                 }
    
        sio.savemat('%(paths.ap)scrab_py' % locals(), p)


###

# EPHEMERIDES
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


###    
        
## SIMULATE            
class Signal(object):
    '''
    Sets absolute properties of a signal: detector, PSR, template and phase difference
    between components. A time vector must also be provided.
    Includes methods to return design matrix and simulated signal given extra inputs
    of polarization and inclination angles.
    '''

    def __init__(self, detector, psr, pdif, t, kind='GR', isloaded=False):
    
        # source and detector information
        self.detector = detector
        self.det = sd.detnames(detector)
        self.psr = psr
        
        # time
        self.t = t
        
        # signal info
        self.kind = kind
        self.basis = sd.aps[kind]
        self.pdif = sd.phase(pdif)
        
        # get antenna patterns
        self.response = Response(psr, detector, t, kind, loadvectors=True)
        # compute roemer delay
        self.getroemer(self, self.response.obs.t)
        
        # parameters
        self.psi = 0
        self.iota = 0
        self.phi0 = 0
        
        self.design_matrix_phase()
            
            
    # Amplitude
    def signalinfo(self, pdif, iota=[], p=0, h_s=0, pdif_s=0):
        
        # determine inclination angle
        if iota==[]:
            # no preset iota, take default
            iota = self.response.src.param['INC']

        # return dataframe with amplitude (h) and phase (phi) for each component
        info = pd.DataFrame(index=self.basis, columns=['h', 'p'])      

        # set amplitudes and phases
        if self.kind == 'GR':
            info['h']['pl'] = (1. + np.cos(iota)**2)/2.
            info['h']['cr'] = np.cos(iota)
            
        elif self.kind == 'G4v':
            info['h']['xz'] = np.sin(iota)
            info['h']['yz'] = np.sin(iota) * math.cos(iota)
            
        elif self.kind == 'GRs':
            info['h']['pl'] = (1. + np.cos(iota)**2)/2.
            info['h']['cr'] = np.cos(iota)
            info['h']['br'] = h_s
                        
        else:
            print 'Warning: %s is not a recognized template (temp 544)'
            exit()
            
        return info
        
        
    def design_matrix_amplitude(self, pol_angle, incl_angle, h_scalar=0):        
        
        if self.kind=='Sid':
            # build basis set:
            th = pd.Series(sd.w * np.array(self.t), index=self.t)
            
            basis = [th.map(np.cos), (2.*th).map(np.cos), th.map(np.sin), (2.*th).map(np.sin)]
            
            # construct matrix
            self.dmA = pd.concat(basis, axis=1, keys=['cos1', 'cos2', 'sin1', 'sin2'])
            self.dmA['cnst']=1

        else:
            print 'creating DMa'
            # make copy of response
            response_local = copy.copy(self.response)
            # get antenna patterns:
            response_local.create(psi=pol_angle, savefile=False)
        
            # get amplitude info
            info = self.signalinfo(0, iota=incl_angle, h_s=h_scalar)
        
            # construct matrix
            dmDict = {pol: getattr(response_local, pol) * info['h'][pol] / 2. for pol in self.basis}
            self.dmA = pd.DataFrame(dmDict).dropna(axis=1) # DF cols: comp. names, index: t.
            self.psi = pol_angle
            self.iota = incl_angle
    
    
    # Phase
    def getroemer(self, t, c=sd.c):
        
        coords = ['x', 'y', 'z']
        
        # unit vector pointing to source (opposite of source-SSB):
        # building following Edwards et al. TEMPO2 (14)
        n0 = - self.response.src.wz
        
        aL = [
                -np.sin(self.response.src.param['RAS']),
                np.cos(self.response.src.param['DEC']),
                0.
            ]
        
        dL = [
                -np.cos(self.response.src.param['RAS'])*np.sin(self.response.src.param['DEC']),
                -np.sin(self.response.src.param['RAS'])*np.sin(self.response.src.param['DEC']),
                np.cos(self.response.src.param['DEC'])
            ]

        a = pd.Series(aL, index=coords)
        d = pd.Series(dL, index=coords)

       
        mu_perp = self.response.src.param['PMRAS']*a + self.response.src.param['PMDEC']*d
        
        dt = pd.Series(self.response.obs.t - self.response.src.param['POSEPOCH'], index=self.response.obs.t)
        
        speed_term = pd.DataFrame(np.outer(mu_perp, dt), index=coords, columns=self.response.obs.t)
        
        accel_term_A = - abs(mu_perp.dot(mu_perp))**2 * n0 /2. 
        accel_term = pd.DataFrame(np.outer(accel_term_A, dt**2), index=coords, columns=self.response.obs.t)
        
        n = (speed_term + accel_term).add(n0, axis='index').T
        
        # vector from geocenter to detector
        r_d = self.response.obs.dz
        
        # vector from SSB to geocenter
        if self.psr == 'J0534+2200':
            ephID = 200
        else:
            ephID = 405
        
        eph = EphemerisD(ephID, self.response.obs.t)
        
        # vector from SSB to detector
        r = r_d + eph.r.T
        
        # dot product
        rn = (r*n).sum(axis=1)
        
        self.roemer = rn/c
    
    def design_matrix_phase(self):

        print 'creating dmP'
        # time
        t = self.response.obs.t
        t2 = t**2
    
        # get frequency and its derivatives
        f = {n : self.response.src.param['F' + str(n)] for n in range(3)}
    
        # get roemer and multiples
        dR = {n : self.roemer ** n for n in range(1,4)}
    
        # form basis vectors
        d3 = dR[3] * f[2] /6.
        d2 = dR[2] * (f[1] + f[2]*t)/2.
        d1 = dR[1] * (f[0] + f[1]*t + f[2]*t**2. /2.)
        d0 = - d1 - d2 - d3
            
        self.dmP = 4. * np.pi * pd.DataFrame({'d3':d3, 'd2':d2, 'd1':d1, 'd0':d0}, index=t)
      
        
    def phase_ev(self, delta, phi0):
        d = pd.Series({'d' + str(n) : delta ** n for n in range(4)})

        dmP = self.dmP.mul(d)               # obtain phase design matrix
        dmP['d0'] += 2. * phi0              # add phi0
        phi_t = dmP.sum(axis='columns')     # sum columns to obtain evolution
        return phi_t
           
    def simulate(self, delta, pol_angle, incl_angle, phase=0, h_scalar=0, pdif_scalar=0, dontadd=False):
        '''
        Simulates a signal based given polarization and inclination angles.
        Includes phase evolution: delta = c/c_g
        Warning: Does not scale output signal, need to multiply output by h0.
        '''
        
        # get amplitude design matrix
        if self.psi==pol_angle and self.iota==incl_angle and 'dmA' in dir(self):
            dm = self.dmA
        else:
            self.design_matrix_amplitude(pol_angle, incl_angle, h_scalar=h_scalar)
            dm = self.dmA

        # construct phase evolution
        phi_t = self.phase_ev(delta, phase)
        
        # compute phase factor
        phase_factor = np.exp(1j*phi_t)
        
        # multiply amplitudes with corresponding phases (mul)
        s1 = dm[self.basis[0]] * phase_factor                            # plus
        
        s2 = dm[self.basis[1]]* phase_factor * np.exp(1j*self.pdif)      # cross
        
        if dontadd:
            return pd.DataFrame({'pl': s1, 'cr': s2})
        else:
            # add columns to form signal
            return s1 + s2