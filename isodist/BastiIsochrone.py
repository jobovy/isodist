import os, os.path
import glob
import csv
import math
import numpy
from Isochrone import Isochrone, FEH2Z, Z2FEH, dict2recarray
from PadovaIsochrone import _DATADIR
_BASTIZSOLAR= 0.0198
_ZS= [0.0001,0.0003]
_ZDICT= {}
_YDICT= {}
_ZDICT['0.0001']= '104'
_YDICT['0.0001']= '245'
_ZDICT['0.0003']= '304'
_YDICT['0.0003']= '245'
class BastiIsochrone (Isochrone):
    """Class that represents a Basti isochrone"""
    def __init__(self,Z=None,filters=None,eta=0.4):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           Z= load only this metallicity (can be list)
           filters= list of filters to load (e.g., ['U','B','V','R','I','J','K','L'])
           eta= (0.4) mass-loss parameter
        OUTPUT:
        HISTORY:
           2012-07-23 - Written - Bovy (IAS)
        """
        if filters is None:
            self._filters= ['U','B','V','R','I','J','H','K','L']
        else:
            self._filters= filters
        #Read the files
        dicts= []
        if Z is None:
            ZS= _ZS
        else:
            if isinstance(Z,(list,numpy.ndarray)):
                ZS= Z
            else:
                ZS= [Z]
        print ZS
        for Zm in ZS:
            if 'K' in self._filters:
                subdir= 'basti-scaled-canonical-%.1f-UBVRIJHKL' % eta
            ages, rawages= _get_ages(os.path.join(_DATADIR,
                                                  subdir,
                                                  'wz'+_ZDICT['%.4f' % Zm]
                                                  +'y'+_YDICT['%.4f' % Zm]
                                                  +'ss2.'+'*'
                                                  +'_c03hbs'))
            dicts.append(read_basti_isochrone(os.path.join(_DATADIR,
                                                           subdir),
                                              'wz'+_ZDICT['%.4f' % Zm]
                                              +'y'+_YDICT['%.4f' % Zm]
                                              +'ss2.',
                                              '_c03hbs',
                                              ages=ages,
                                              rawages=rawages,
                                              filters=self._filters))
        self._ZS= ZS
        self._dicts= dicts
        #Gather ages
        self._logages= numpy.array(sorted(list(set(self._dicts[0]['logage']))))
        return None
        
def _get_ages(path):
    """Return the available ages for this isochrone set"""
    files= glob.glob(path)
    rawages= []
    for ii in range(len(files)):
        rawages.append((((os.path.basename(files[ii]).split('.'))[1]).split('_'))[0])
    return (numpy.array([_parse_age(r) for r in rawages]),rawages)

def _parse_age(raw):
    """Parse a raw age (e.g., 't600030')"""
    return float(raw[2:])*10.**-3. #In Gyr

def read_basti_isochrone(dir,name1,name2,ages=None,rawages=None,
                         filters=None):
    """
    NAME:
       read_basti_isochrone
    PURPOSE:
       read a Basti isochrone file
    INPUT:
       dir - directory of the filename
       name1 - first part of the name (before the age)
       name2 - second part of the name(after the age)
       filters= list of filters in the file
       age= age in Gyr
       rawages= raw age strings for the filenames
    OUTPUT:
       dictionary with the table
    HISTORY:
       2012-07-23 - Written - Bovy (IAS)
    """
    nages= len(ages)
    dialect= csv.excel
    dialect.skipinitialspace=True
    nfilters= len(filters)
    ncols= nfilters+8
    logage=[]
    y= []
    M_ini= []
    M_act= []
    logL= []
    logTe= []
    mags= []
    for ii in range(nages):
        file= open(os.path.join(dir,name1+rawages[ii]+name2),'r')
        reader= csv.reader(file,delimiter=' ',
                           dialect=dialect)
        for row in reader:
            try:
                if row[0][0] == '#':
                    continue
            except IndexError:
                pass
            logage.append(9.+math.log10(ages[ii]))
            M_ini.append(float(row[0]))
            M_act.append(float(row[1]))
            logL.append(float(row[2]))
            logTe.append(float(row[3]))
            if 'K' in filters:
                V= float(row[4])
                B= float(row[6])+V
                U= float(row[5])+B
                I= V-float(row[7])
                R= V-float(row[8])
                J= V-float(row[9])
                K= V-float(row[10])
                L= V-float(row[11])
                H= float(row[12])+K
            mags.append([U,B,V,R,I,J,H,K,L])
    #Load everything into a dictionary
    outDict= {}
    outDict['logage']= numpy.array(logage)
    outDict['M_ini']= numpy.array(M_ini)
    outDict['M_act']= numpy.array(M_act)
    outDict['logL']= numpy.array(logL)
    outDict['logTe']= numpy.array(logTe)
    for ii in range(nfilters):
        thismag= []
        for jj in range(len(mags)):
            thismag.append(mags[jj][ii])
        outDict[filters[ii]]= numpy.array(thismag)
    return outDict
