import os, os.path
import glob
import csv
import math
import numpy
from Isochrone import Isochrone, FEH2Z, Z2FEH, dict2recarray, logg
from PadovaIsochrone import _DATADIR
_BASTIZSOLAR= 0.0198
_ZS= [0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,0.01,0.0198,
      0.03,0.04]
_ZDICT= {}
_YDICT= {}
_ZDICT['0.0001']= '104'
_YDICT['0.0001']= '245'
_ZDICT['0.0003']= '304'
_YDICT['0.0003']= '245'
_ZDICT['0.0006']= '604'
_YDICT['0.0006']= '246'
_ZDICT['0.0010']= '103'
_YDICT['0.0010']= '246'
_ZDICT['0.0020']= '203'
_YDICT['0.0020']= '248'
_ZDICT['0.0040']= '403'
_YDICT['0.0040']= '251'
_ZDICT['0.0080']= '803'
_YDICT['0.0080']= '256'
_ZDICT['0.0100']= '102'
_YDICT['0.0100']= '259'
_ZDICT['0.0198']= 'sun'
_YDICT['0.0198']= 'sun'
_ZDICT['0.0300']= '302'
_YDICT['0.0300']= '288'
_ZDICT['0.0400']= '402'
_YDICT['0.0400']= '303'
#Dictionary for last part of filename
post= {}
post['UBVRIJHKL']= 'c03hbs'
post['ugriz']= 'sloan'
post['uu_0bym1c1c1_0H_betaCa']= 'strm'
strmfilters= ['u','u_0','b','y','m1','c1','c1_0','H_beta','Ca']
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
        for Zm in ZS:
            subdir= 'basti-scaled-canonical-%.1f-' % eta + ''.join(self._filters)
            if ''.join(self._filters) == 'uu_0bym1c1c1_0H_betaCa':
                subdir= 'basti-scaled-canonical-%.1f-strm' % eta
            if eta == 0.2:
                etastr= 'ss2.'
            elif eta == 0.4:
                etastr= 's.'
            ages, rawages= _get_ages(os.path.join(_DATADIR,
                                                  subdir,
                                                  'wz'+_ZDICT['%.4f' % Zm]
                                                  +'y'+_YDICT['%.4f' % Zm]
                                                  +etastr+'*'
                                                  +'_'+post[''.join(self._filters)]))
            dicts.append(read_basti_isochrone(os.path.join(_DATADIR,
                                                           subdir),
                                              'wz'+_ZDICT['%.4f' % Zm]
                                              +'y'+_YDICT['%.4f' % Zm]
                                              +etastr,
                                              '_'+post[''.join(self._filters)],
                                              ages=ages,
                                              rawages=rawages,
                                              filters=self._filters))
        self._ZS= numpy.array(ZS)
        self._dicts= dicts
        #Gather ages
        self._logages= numpy.array(sorted(list(set(self._dicts[0]['logage']))))
        return None
        
    def __call__(self,logage,Z=None,feh=None,afe=None,maxm=None,
                 asrecarray=False,stage=None):
        """
        NAME:
           __call__
        PURPOSE:
           get a single isochrone from the library
        INPUT:
           logage - log_10 age
           Z= or feh= metallicity (use Z_\odot=0.019)
           afe= None (not yet supported for Basti)
           maxm= maximum mass to consider (m_ini)
           stage= if set, only show this evolutionary stage (NOT IMPLEMENTED FOR BASTI)
        KEYWORDS:
           asrecarray= if True, return recarray
        OUTPUT:
           isochrone
        HISTORY:
           2012-07-23 - Written - Bovy (IAS)
        """
        if not afe is None:
            raise NotImplementedError("'afe=' not yet implemented for Basti isochrones")
        if not feh is None:
            Z= FEH2Z(feh)
        indx= (self._ZS == Z)
        ii= 0
        while (not indx[ii]): ii+= 1
        thisDict= self._dicts[ii]
        if maxm is None:
            indx= (thisDict['logage'] == logage)
        else:
            indx= (thisDict['logage'] == logage)*(thisDict['M_ini'] < maxm)
        outDict= {}
        for key in thisDict.keys():
            outDict[key]= thisDict[key][indx]
        if asrecarray:
            return dict2recarray(outDict)
        else:
            return outDict

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
            if ''.join(filters) == 'UBVRIJHKL':
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
            elif ''.join(filters) == 'ugriz':
                g= float(row[4])
                u= float(row[5])+g
                r= g-float(row[6])
                i= r-float(row[7])
                z= i-float(row[8])
                mags.append([u,g,r,i,z])
            elif ''.join(filters) == 'uu_0bym1c1c1_0H_betaCa':
                y= float(row[4])
                b= y+float(row[7])
                u= b+float(row[5])
                u_0= b+float(row[6])
                m1= float(row[8])
                c1= float(row[9])
                c1_0= float(row[10])
                H_beta= float(row[11])
                Ca= float(row[12])
                mags.append([u,u_0,b,y,m1,c1,c1_0,H_beta,Ca])
    #Load everything into a dictionary
    outDict= {}
    outDict['logage']= numpy.array(logage)
    outDict['M_ini']= numpy.array(M_ini)
    outDict['M_act']= numpy.array(M_act)
    outDict['logL']= numpy.array(logL)
    outDict['logTe']= numpy.array(logTe)
    outDict['logg']= logg(outDict['logL'],outDict['logTe'],outDict['M_act'])
    for ii in range(nfilters):
        thismag= []
        for jj in range(len(mags)):
            thismag.append(mags[jj][ii])
        outDict[filters[ii]]= numpy.array(thismag)
    return outDict
