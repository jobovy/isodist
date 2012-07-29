###############################################################################
#   DartmouthIsochrone: Module that represents isochrones calculated by the 
#                       Dartmouth group
#
#   Quick start guide (NEEDS TO BE EDITED FOR DARTMOUTH)
#   -----------------
#
#   - Initialization: >>>p= PadovaIsochrone()
#                     or
#                     >>>p= PadovaIsochrone(Z=[0.02,0.03],
#                                           type='2mass-spitzer-wise')
#
#   - Calling: >>>p(8.,Z=0.02)
#              returns the whole isochrone as a dictionary
#
#   - Plotting: >>>p.plot(logage,Z=,feh=,d1=,d2=,maxm=,...)
#               E.g.,
#               >>>p.plot(8.,Z=0.02,d1='J-Ks',d2='H')
#               plots the log_10 age= 8., Z=0.02 isochrone in the J-Ks,H plane
#
#   - Defining a new PadovaIsochrone: 1) In $ISODIST_DATA make a new folder 
#                                        with the name you want (e.g., '2mass')
#                                     2) put lists of isochrones of a certain Z
#                                        and age (one file / Z) in this 
#                                        directory; these files are named 
#                                        __NAME__-Z-__Z__+.gz
#                                        e.g.,
#                                        2mass-Z-0.020.dat.gz
#                                     3) Edit the __init__ of the 
#                                        PadovaIsochrone to include this new 
#                                        type (optional)
###############################################################################
import os, os.path
import csv
import gzip
import math
import numpy
from Isochrone import Isochrone, FEH2Z, Z2FEH, dict2recarray
_FEHS= [-2.5,-2.,-1.5,-1.,-0.5,0.,0.2,0.3,0.5]
_DATADIR= os.getenv('ISODIST_DATA')
#Dictionary for last part of filename
post= {}
post['UBVRIJHKs']= 'jc2mass'
if _DATADIR is None:
    _DATADIR= os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           '../data')
class DartmouthIsochrone (Isochrone):
    """Class that represents a Dartmouth isochrone"""
    def __init__(self,feh=None,filters=None,afe=0.):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           Z= load only this metallicity (can be list)
           filters= list of filters (optional)
           afe= [a/Fe] (default: 0.)
        OUTPUT:
        HISTORY:
           2012-07-29 - Written - Bovy (IAS)
        """
        #Set filters
        if filters is None:
            self._filters= ['U','B','V','R','I','J','H','Ks']
        else:
            self._filters= filters
        #Read the files
        dicts= []
        if feh is None:
            FEHS= _FEHS
        else:
            if isinstance(feh,(list,numpy.ndarray)):
                FEHS= feh
            else:
                FEHS= [feh]
        for fehm in FEHS:
            if fehm >= 0.: fehsignstr= 'p'
            else: fehsignstr= 'm'
            if afe >= 0.: afesignstr= 'p'
            else: afesignstr= 'm'
            dicts.append(read_dartmouth_isochrone(os.path.join(_DATADIR,
                                                               'dartmouth-'+"".join(self._filters),
                                                               'feh'+fehsignstr+'%02i' % (int(numpy.fabs(10.*fehm)))\
                                                                   +'afe'+afesignstr+'%01i' % (int(numpy.fabs(10.*afe)))\
                                                                   +'.'+post["".join(self._filters)]),
                                                               filters=self._filters))
        self._ZS= FEH2Z(numpy.array(FEHS))
        self._dicts= dicts
        #Gather ages
        self._logages= numpy.array(sorted(list(set(self._dicts[0]['logage']))))
        return None

    def __call__(self,logage,Z=None,feh=None,afe=None,maxm=None,
                 asrecarray=False):
        """
        NAME:
           __call__
        PURPOSE:
           get a single isochrone from the library
        INPUT:
           logage - log_10 age
           Z= or feh= metallicity (use Z_\odot=0.019)
           afe= None (not supported for Padova)
           maxm= maximum mass to consider (m_ini)
        KEYWORDS:
           asrecarray= if True, return recarray
        OUTPUT:
           isochrone
        HISTORY:
           Written - Bovy (NYU)
        """
        if not afe is None:
            raise NotImplementedError("'afe=' not implemented for Padova isochrones")
        if not feh is None:
            Z= FEH2Z(feh)
        indx= (self._ZS == Z)
        ii= 0
        while (not indx[ii]): ii+= 1
        thisDict= self._dicts[ii]
        if maxm is None:
            indx= (thisDict['logage'] == logage)
        else:
            indx= (thisDict['logage'] == logage)*(thisDict['M'] < maxm)
        outDict= {}
        for key in thisDict.keys():
            outDict[key]= thisDict[key][indx]
        if asrecarray:
            return dict2recarray(outDict)
        else:
            return outDict

def read_dartmouth_isochrone(name,filters=None,skiponegyr=False):
    """
    NAME:
       read_dartmouth_isochrone
    PURPOSE:
       read a Dartmouth isochrone file
    INPUT:
       name- name of the file
       filters= list of filters in the file
       skiponegyr= if True, skip the 1 Gyr isochrone (since it's both in the young and old isochrones)
    OUTPUT:
       dictionary with the table
    HISTORY:
       2012-07-29 - Written - Bovy (IAS)
    """
    dialect= csv.excel
    dialect.skipinitialspace=True
    files= [open(name,'r'),open(name+'_2','r')]
    nfilters= len(filters)
    ncols= nfilters+5
    logage=[]
    EEP= []
    M= []
    logL= []
    logTe= []
    logg= []
    mags= []
    second= False
    for file in files:
        reader= csv.reader(file,delimiter=' ',
                           dialect=dialect)
        for row in reader:
            try:
                if row[0][0] == '#':
                    if row[0][1:4] == 'AGE':
                        try:
                            currentage= float(row[1])
                        except ValueError:
                            currentage= float(row[0].split('=')[1])
                    continue
            except IndexError:
                if len(row) == 0: continue
                pass
            if second and currentage == 1.: continue
            logage.append(6.+numpy.log10(currentage))
            EEP.append(float(row[0]))
            M.append(float(row[1]))
            logTe.append(float(row[2]))
            logg.append(float(row[3]))
            logL.append(float(row[4]))
            thismags= []
            for ii in range(nfilters):
                thismags.append(float(row[5+ii]))
            mags.append(thismags)
            second= True
    #Load everything into a dictionary
    outDict= {}
    outDict['logage']= numpy.array(logage)
    outDict['EEP']= numpy.array(EEP)
    outDict['M']= numpy.array(M)
    outDict['logL']= numpy.array(logL)
    outDict['logTe']= numpy.array(logTe)
    outDict['logg']= numpy.array(logg)
    for ii in range(nfilters):
        thismag= []
        for jj in range(len(mags)):
            thismag.append(mags[jj][ii])
        outDict[filters[ii]]= numpy.array(thismag)
    return outDict
