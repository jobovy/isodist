import os, os.path
import csv
import gzip
import numpy as nu
_ZS= [0.002,0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02,0.022,
      0.024,0.026,0.028,0.03]
_DATADIR= os.getenv('ISODIST_DATA')
if _DATADIR is None:
    _DATADIR= '../data'
class PadovaIsochrone:
    """Class that represents a Padova isochrone"""
    def __init__(self,type='2mass-spitzer-wise',Z=None):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           type= type of isochrones to load (e.g., 2mass-spitzer-wise)
           Z= load only this metallicity (can be list)
        OUTPUT:
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        #Set filters
        if type.lower() == '2mass-spitzer-wise':
            filters= ['J','H','Ks','[3.6]','[4.5]','[5.8]','[8.0]','[24]','[70]','[160]','W1','W2','W3','W4']
        #Read the files
        dicts= []
        if Z is None:
            ZS= _ZS
        else:
            if isinstance(Z,list):
                ZS= Z
            else:
                ZS= [Z]
        for Zm in ZS:
            dicts.append(read_padova_isochrone(os.path.join(_DATADIR,
                                                            type.lower(),
                                                            type.lower()+'-Z-%5.3f.dat.gz' % Zm),
                                               filters=filters))
        self._dicts= dicts
        return None

def read_padova_isochrone(name,filters=None):
    """
    NAME:
       read_padova_isochrone
    PURPOSE:
       read a Padova isochrone file
    INPUT:
       name- name of the file
       filters= list of filters in the file
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
        file= gzip.open(name,'r')
    else:
        file= open(name,'r')
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
    CO= []
    M_hec= []
    period= []
    pmode= []
    logMdot= []
    int_IMF= []
    mags= []
    for row in reader:
        try:
            if row[0][0] == '#':
                continue
        except IndexError:
            pass
        logage.append(float(row[1]))
        M_ini.append(float(row[2]))
        M_act.append(float(row[3]))
        logL.append(float(row[4]))
        logTe.append(float(row[5]))
        logg.append(float(row[6]))
        mbol.append(float(row[7]))
        thismags= []
        for ii in range(nfilters):
            thismags.append(float(row[8+ii]))
        mags.append(thismags)
        CO.append(float(row[8+nfilters]))
        M_hec.append(float(row[9+nfilters]))
        period.append(float(row[10+nfilters]))
        pmode.append(float(row[11+nfilters]))
        logMdot.append(float(row[12+nfilters]))
        int_IMF.append(float(row[13+nfilters]))
    #Load everything into a dictionary
    outDict= {}
    outDict['logage']= nu.array(logage)
    outDict['M_ini']= nu.array(M_ini)
    outDict['M_act']= nu.array(M_act)
    outDict['logL']= nu.array(logL)
    outDict['logTe']= nu.array('logTe')
    outDict['logg']= nu.array(logg)
    outDict['mbol']= nu.array(mbol)
    outDict['CO']= nu.array(CO)
    outDict['M_hec']= nu.array(M_hec)
    outDict['period']= nu.array(period)
    outDict['pmode']= nu.array(pmode)
    outDict['logMdot']= nu.array(logMdot)
    outDict['int_IMF']= nu.array(int_IMF)
    for ii in range(nfilters):
        thismag= []
        for jj in range(len(mags)):
            thismag.append(mags[jj][ii])
        outDict[filters[ii]]= nu.array(thismag)
    return outDict
