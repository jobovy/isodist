import os, os.path
import re
import csv
import gzip
import math
import numpy as nu
from galpy.util import bovy_plot
_ZS= [0.002,0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02,0.022,
      0.024,0.026,0.028,0.03]
_DATADIR= os.getenv('ISODIST_DATA')
_ZSOLAR= 0.019
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
        self._ZS= nu.array(ZS)
        self._dicts= dicts
        self._filters= filters
        return None

    def __call__(self,logage,Z=None,feh=None,afe=None,maxm=None):
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
        OUTPUT:
           isochrone
        HISTORY:
           Written - Bovy (NYU)
        """
        if not afe is None:
            raise NotImplementedError("'afe=' not implemented for Padova isochrones")
        if not feh is None:
            Z= math.exp(feh+math.log(_ZSOLAR))
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
        outDict= {}
        for key in thisDict.keys():
            outDict[key]= thisDict[key][indx]
        return outDict

    def plot(self,logage,*args,**kwargs):
        """
        NAME:
           plot
        PURPOSE:
           plot an individual isochrone, or a set of isochrones
        INPUT:
           logage - log_10 age or list thereof
           Z= metallicity or list
           feh= metallicity or list (use Z_\odot=0.019)
           afe= not supported for Padova
           d1= x dimension (for color write 'J-Ks')
           d2= y dimension
           +bovy_plot.bovy_plot keywords
        OUTPUT:
           plot to output device
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        if not isinstance(logage,(list,nu.ndarray)):
            return self._plot_single(logage,*args,**kwargs)
        #Handle Z etc.
        if kwargs.has_key('feh'):
            if isinstance(kwargs['feh'],(list,nu.ndarray)):
                fehs= kwargs['feh']
            else:
                fehs= [kwargs['feh'] for ii in range(len(logage))]
        else:
            fehs= [None for ii in range(len(logage))]
        if kwargs.has_key('Z'):
            if isinstance(kwargs['Z'],(list,nu.ndarray)):
                ZS= kwargs['Z']
            else:
                ZS= [kwargs['Z'] for ii in range(len(logage))]
        else: ZS= [None for ii in range(len(logage))]
        #Plot first
        out= self._plot_single(logage[0],*args,**kwargs)
        if not kwargs.has_key('overplot') \
                or (kwargs.has_key('overplot') and not kwargs['overplot']):
            kwargs['overplot']= True
        for ii in range(1,len(logage)):
            kwargs['Z']= ZS[ii]
            kwargs['feh']= fehs[ii]
            self._plot_single(logage[ii],*args,**kwargs)
        return out

    def _plot_single(self,logage,*args,**kwargs):
        #kwargs
        if kwargs.has_key('Z'):
            Z= kwargs['Z']
            kwargs.pop('Z')
        else: Z= None
        if kwargs.has_key('feh'):
            feh= kwargs['feh']
            kwargs.pop('feh')
        else: feh= None
        if kwargs.has_key('afe'):
            afe= kwargs['afe']
            kwargs.pop('afe')
        else: afe= None
        if kwargs.has_key('maxm'):
            maxm= kwargs['maxm']
            kwargs.pop('maxm')
        else: maxm= None
        if kwargs.has_key('d1'):
            d1= kwargs['d1']
            kwargs.pop('d1')
        else: d1= self._filters[0]+'-'+self._filters[1]
        if kwargs.has_key('d2'):
            d2= kwargs['d2']
            kwargs.pop('d2')
        else: d2= self._filters[0]
        #get isochrone
        iso= self(logage,Z=Z,feh=feh,afe=afe,maxm=maxm)
        #get dimensions
        colorx= re.split(r'-',d1)
        if len(colorx) == 2: #d1 is a color
            x= iso[colorx[0]]-iso[colorx[1]]
        else:
            x= iso[d1]
        colory= re.split(r'-',d2)
        if len(colory) == 2: #d2 is a color
            y= iso[colory[0]]-iso[colory[1]]
        else:
            y= iso[d2]
        #Put in default labels
        if not kwargs.has_key('overplot') \
                or (kwargs.has_key('overplot') and not kwargs['overplot']):
            if not kwargs.has_key('xlabel'):
                kwargs['xlabel']= r'$'+d1+'$'
            if not kwargs.has_key('ylabel'):
                if d2 in self._filters:
                    kwargs['ylabel']= r'$M_{'+d2+'}$'
                else:
                    kwargs['ylabel']= r'$'+d2+'$'
            if not kwargs.has_key('yrange') and d2 in self._filters:
                kwargs['yrange']= [nu.amax(y)+0.3,nu.amin(y)-0.3]
        #plot
        return bovy_plot.bovy_plot(x,y,*args,**kwargs)

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
    outDict['logTe']= nu.array(logTe)
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
