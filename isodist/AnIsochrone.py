import os, os.path
import csv
import math
import numpy
from Isochrone import Isochrone, FEH2Z, Z2FEH, dict2recarray
from PadovaIsochrone import _DATADIR
_ANZSOLAR= 0.0176
_ZS= [-0.1,-0.2,-0.3,-0.5,-1.,-1.5,-2.,-3.,0.,0.1,0.2,0.4]
class AnIsochrone (Isochrone):
    """Class that represents a An+08 isochrone"""
    def __init__(self,Z=None,filters=None,corrected=True):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           corrected= if False, use un-corrected isochrones
           Z= load only this metallicity (can be list)
        OUTPUT:
        HISTORY:
           2011-08-05 - Written - Bovy (NYU)
        BUGS:
           Z determination needs to account for dY/dZ
        """
        self._filters= ['u','g','r','i','z']
        #Read the files
        dicts= []
        if Z is None: #Z here is actually FeH, we correct this later
            ZS= _ZS
        else:
            if isinstance(Z,(list,numpy.ndarray)):
                ZS= Z
            else:
                ZS= [Z]
        for Zm in ZS:
            if Zm >= 0.: signstr= 'p'
            else: signstr= 'm'
            if corrected: corrstr= 'corr'
            else: corrstr= 'marcs'
            dicts.append(read_an_isochrone(os.path.join(_DATADIR,
                                                        'an_isochrones',
                                                        signstr+'%03i_' % (int(numpy.fabs(100.*Zm)))
                                                        +corrstr+'.txt'),
                                           filters=self._filters))
        self._ZS= numpy.array([FEH2Z(z,zsolar=_ANZSOLAR) for z in ZS])
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
           afe= None (not supported for An; linear relation between afe and 
                feh is assumed)
           maxm= maximum mass to consider (m_ini)
        KEYWORDS:
           asrecarray= if True, return recarray, otherwise dict
        OUTPUT:
           isochrone
        HISTORY:
           2011-08-04 - Written - Bovy (NYU)
        """
        if not afe is None:
            raise NotImplementedError("'afe=' not implemented for Padova isochrones")
        if not feh is None:
            Z= math.exp(feh+math.log(_ANZSOLAR))
        indx= (self._ZS == Z)
        ii= 0
        while (ii < len(self._dicts) and not indx[ii]): ii+= 1
        if ii == len(self._dicts):
            raise IOError("No isochrone found that matches this metallicity")
        thisDict= self._dicts[ii]
        if maxm is None:
            indx= (thisDict['logage'] == logage)
        else:
            indx= (thisDict['logage'] == logage)*(thisDict['Mass'] < maxm)
        if numpy.sum(indx) == 0:
            raise IOError("No isochrone found that matches this logage")
        outDict= {}
        for key in thisDict.keys():
            outDict[key]= thisDict[key][indx]
        if asrecarray:
            return dict2recarray(outDict)
        else:
            return outDict

def read_an_isochrone(name,filters=None):
    """
    NAME:
       read_an_isochrone
    PURPOSE:
       read an An isochrone file
    INPUT:
       name- name of the file
       filters= list of filters in the file
    OUTPUT:
       dictionary with the table
    HISTORY:
       2011-08-04 - Written - Bovy (NYU)
    """
    dialect= csv.excel
    dialect.skipinitialspace=True
    if name[-2:] == 'gz':
        file= gzip.open(name,'r')
    else:
        file= open(name,'r')
    reader= csv.reader(file,delimiter=' ',
                       dialect=dialect)
    ncols= 10
    logage=[]
    Mass= []
    logL= []
    logTe= []
    logg= []
    mbol= []
    mags= []
    for row in reader:
        try:
            if row[0][0:4] == 'Mass': #Header line to skip
                continue
        except IndexError:
            pass
        try:
            if row[0] == 'Cluster': #Header line to extract age from
                thislogage= numpy.log10(float(row[4]))
                continue
        except IndexError:
            pass
        logage.append(thislogage) #from the header, see above
        Mass.append(float(row[0]))
        logTe.append(numpy.log10(float(row[1])))
        logL.append(float(row[2]))
        logg.append(float(row[3]))
        mbol.append(float(row[4]))
        r= float(row[5])
        gr = float(row[6])
        gi = float(row[7])
        gz = float(row[8])
        ug = float(row[9])
        mags.append([r+gr+ug, #u
                     r+gr, #g
                     r,
                     -gi+gr+r, #i
                     -gz+gr+r]) #z
    #Load everything into a dictionary
    outDict= {}
    outDict['logage']= numpy.array(logage)
    outDict['Mass']= numpy.array(Mass)
    outDict['logL']= numpy.array(logL)
    outDict['logTe']= numpy.array(logTe)
    outDict['logg']= numpy.array(logg)
    outDict['mbol']= numpy.array(mbol)
    for ii in range(len(filters)):
        thismag= []
        for jj in range(len(mags)):
            thismag.append(mags[jj][ii])
        outDict[filters[ii]]= numpy.array(thismag)
    return outDict
