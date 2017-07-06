###############################################################################
#   PadovaIsochrone: Module that represents isochrones calculated by the 
#                    Padova group
#
#   Quick start guide
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
import sys
import os, os.path
import csv
import gzip
import math
import numpy as nu
from isodist.Isochrone import Isochrone, FEH2Z, Z2FEH, dict2recarray
_ZS= [0.002,0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02,0.022,
      0.024,0.026,0.028,0.03]
_DATADIR= os.getenv('ISODIST_DATA')
if _DATADIR is None:
    _DATADIR= os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           '../data')
if sys.version_info[0] < 3:
    FileNotFoundError= IOError
class PadovaIsochrone (Isochrone):
    """Class that represents a Padova isochrone"""
    def __init__(self,type='2mass-spitzer-wise',Z=None,filters=None,
                 parsec=False,eta=None):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           type= type of isochrones to load (e.g., 2mass-spitzer-wise)
           Z= load only this metallicity (can be list)
           filters= list of filters (optional)
           parsec= if True, use new PARSEC isochrones
           eta= Reimers mass loss efficiency parameter 
                (default: 0.4 for Padova, 0.2 for PARSEC)
        OUTPUT:
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        #Set filters
        if type.lower() == '2mass-spitzer-wise':
            filters= ['J','H','Ks','[3.6]','[4.5]','[5.8]','[8.0]','[24]','[70]','[160]','W1','W2','W3','W4']
        elif type.lower() == 'sdss-ukidss':
            filters= ['u','g','r','i','z','Z','Y','J','H','K']
        elif type.lower() == 'sdss-2mass':
            filters= ['u','g','r','i','z','J','H','Ks']
        #Read the files
        dicts= []
        if Z is None:
            ZS= _ZS
        else:
            if isinstance(Z,(list,nu.ndarray)):
                ZS= Z
            else:
                ZS= [Z]
        self.parsec= parsec
        if parsec:
            if eta == 0.2 or eta is None:
                basename= 'parsec-'+type.lower()
            else:
                basename= 'parsec-%.1f-' % eta +type.lower()
        else:
            if eta == 0.4 or eta is None:
                basename= type.lower()
            else:
                raise NotImplementedError('Non-default eta not implemented yet for Padova isochrones')
        for Zm in ZS:
            try:
                dicts.append(read_padova_isochrone(os.path.join(_DATADIR,
                                                                basename,
                                                                basename+'-Z-%5.3f.dat.gz' % Zm),
                                                   filters=filters,
                                                   parsec=parsec))
            except FileNotFoundError: #4?
                dicts.append(read_padova_isochrone(os.path.join(_DATADIR,
                                                                basename,
                                                                basename+'-Z-%5.4f.dat.gz' % Zm),
                                                   filters=filters,
                                                   parsec=parsec))
                
        self._ZS= nu.array(ZS)
        self._dicts= dicts
        self._filters= filters
        #Gather ages
        self._logages= nu.array(sorted(list(set(self._dicts[0]['logage']))))
        return None

    def __call__(self,logage,Z=None,feh=None,afe=None,maxm=None,
                 asrecarray=False,
                 stage=None):
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
           stage= if set, only show this evolutionary stage
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
            Z= FEH2Z(feh,parsec=self.parsec)
        indx= (self._ZS == Z)
        ii= 0
        while (not indx[ii]): ii+= 1
        thisDict= self._dicts[ii]
        #round logage
        logage= round(100.*logage)/100.
        if maxm is None:
            indx= (thisDict['logage'] == logage)
        else:
            indx= (thisDict['logage'] == logage)*(thisDict['M_ini'] < maxm)
        if not stage is None:
            indx *= (thisDict['stage'] == stage)
        outDict= {}
        for key in thisDict.keys():
            outDict[key]= thisDict[key][indx]
        if asrecarray:
            return dict2recarray(outDict)
        else:
            return outDict

def read_padova_isochrone(name,filters=None,parsec=False):
    """
    NAME:
       read_padova_isochrone
    PURPOSE:
       read a Padova isochrone file
    INPUT:
       name- name of the file
       filters= list of filters in the file
       parsec= if True, use new parsec isochrones, which have stage
    OUTPUT:
       dictionary with the table
    EXAMPLE:
       read_padova_isochrone('../data/output411373137337.dat',filters=['J','H','Ks','[3.6]','[4.5]','[5.8]','[8.0]','[24]','[70]','[160]','W1','W2','W3','W4'])
    HISTORY:
       2011-04-26 - Written - Bovy (NYU)
    """
    dialect= csv.excel
    dialect.skipinitialspace=True
    if name[-2:] == 'gz':
        file= gzip.open(name,'rt')
    else:
        file= open(name,'rt')
    reader= csv.reader(file,delimiter='\t',
                       dialect=dialect)
    nfilters= len(filters)
    ncols= nfilters+13
    logage=[]
    M_ini= []
    M_act= []
    logL= []
    logTe= []
    logg= []
    mbol= []
    if not parsec:
        CO= []
        M_hec= []
        period= []
        pmode= []
        logMdot= []
    int_IMF= []
    mags= []
    if parsec: stage= []
    for row in reader:
        try:
            if row[0][0] == '#':
                continue
        except IndexError:
            pass
        logage.append(float(row[1+parsec])) #parsec has Z as a column
        M_ini.append(float(row[2+parsec]))
        M_act.append(float(row[3+parsec]))
        logL.append(float(row[4+parsec]))
        logTe.append(float(row[5+parsec]))
        logg.append(float(row[6+parsec]))
        mbol.append(float(row[7+parsec]))
        thismags= []
        for ii in range(nfilters):
            thismags.append(float(row[8+parsec+ii]))
        mags.append(thismags)
        if not parsec:
            CO.append(float(row[8+nfilters]))
            M_hec.append(float(row[9+nfilters]))
            period.append(float(row[10+nfilters]))
            pmode.append(float(row[11+nfilters]))
            logMdot.append(float(row[12+nfilters]))
        int_IMF.append(float(row[13-4*parsec+nfilters])) #4 bc of extra Z
        if parsec:
            stage.append(float(row[14-4*parsec+nfilters]))
    #Load everything into a dictionary
    outDict= {}
    outDict['logage']= nu.array(logage)
    outDict['M_ini']= nu.array(M_ini)
    outDict['M_act']= nu.array(M_act)
    outDict['logL']= nu.array(logL)
    outDict['logTe']= nu.array(logTe)
    outDict['logg']= nu.array(logg)
    outDict['mbol']= nu.array(mbol)
    if not parsec:
        outDict['CO']= nu.array(CO)
        outDict['M_hec']= nu.array(M_hec)
        outDict['period']= nu.array(period)
        outDict['pmode']= nu.array(pmode)
        outDict['logMdot']= nu.array(logMdot)
    outDict['int_IMF']= nu.array(int_IMF)
    if parsec:
        outDict['stage']= nu.array(stage)
    for ii in range(nfilters):
        thismag= []
        for jj in range(len(mags)):
            thismag.append(mags[jj][ii])
        outDict[filters[ii]]= nu.array(thismag)
    return outDict

def padovaTypes():
    return ['2mass-spitzer-wise','sdss-ukidss']
