import re
import math
import numpy
try:
    from galpy.util import bovy_plot
except ImportError:
    import bovy_plot #BOVY: COPY IN LOCAL VERSION
_ZSOLAR= 0.019
_LOGTESUN= numpy.log10(5777)
_LOGGSUN= numpy.log10(27400.)
class Isochrone:
    """Template for any Isochrone type class"""
    def __init__(self):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
        OUTPUT:
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        raise NotImplementedError("'__init__' not implemented for this isochrone")

    def __call__(self,logage,Z=None,feh=None,afe=None,maxm=None,
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
           stage= if set, only return this evolutionary stage 
                  (if this exists for this isochrone libary)
        OUTPUT:
           isochrone
        HISTORY:
           Written - Bovy (NYU)
        """
        raise NotImplementedError("'__call__' method for this isochrone is not implemented")

    def Zs(self):
        """
        NAME:
           Zs
        PURPOSE:
           Return the loaded metallicities
        INPUT:
        OUTPUT:
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        return self._ZS

    def logages(self):
        """
        NAME:
           logages
        PURPOSE:
           Return the loaded log_10 ages
        INPUT:
        OUTPUT:
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        return self._logages

    def filters(self):
        """
        NAME:
           filters
        PURPOSE:
           Return the supported filters
        INPUT:
        OUTPUT:
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        return self._filters
###################################PLOTTING####################################
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
           maxm= maximum mass to plot
           stage= if set, only return this evolutionary stage 
                  (if this exists for this isochrone libary)
           ignore_gaps= if True, ignore non-existant isochrones
           +bovy_plot.bovy_plot keywords
        OUTPUT:
           plot to output device
        HISTORY:
           2011-04-27 - Written - Bovy (NYU)
        """
        if not isinstance(logage,(list,numpy.ndarray)) \
                and not ((kwargs.has_key('Z') \
                              and isinstance(kwargs['Z'],(list,numpy.ndarray)))\
                             or (kwargs.has_key('feh') \
                                     and isinstance(kwargs['feh'],numpy.ndarray))):
            return self._plot_single(logage,*args,**kwargs)
        #Do we have Z or FeH?
        if not kwargs.has_key('Z') and kwargs.has_key('feh'): usefeh= True
        else: usefeh= False
        if not isinstance(logage,(list,numpy.ndarray)) and usefeh:
            logage= numpy.array([logage for ii in range(len(kwargs['feh']))])
        elif not isinstance(logage,(list,numpy.ndarray)):
            logage= numpy.array([logage for ii in range(len(kwargs['Z']))])
        #Handle Z etc.
        if kwargs.has_key('feh'):
            if isinstance(kwargs['feh'],(list,numpy.ndarray)):
                fehs= kwargs['feh']
            else:
                fehs= numpy.array([kwargs['feh'] for ii in range(len(logage))])
        else:
            fehs= [None for ii in range(len(logage))]
        if kwargs.has_key('Z'):
            if isinstance(kwargs['Z'],(list,numpy.ndarray)):
                ZS= kwargs['Z']
            else:
                ZS= numpy.array([kwargs['Z'] for ii in range(len(logage))])
        else: ZS= [None for ii in range(len(logage))]
        if usefeh and not len(logage) == len(fehs):
            raise IOError("When both logage and feh are given as arrays they need to have the same length")
        elif not usefeh and not len(logage) == len(ZS):
            raise IOError("When both logage and Z are given as arrays they need to have the same length")
        #Plot first
        if usefeh: kwargs['feh']= fehs[0]
        else: kwargs['Z']= ZS[0]
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
        if kwargs.has_key('stage'):
            stage= kwargs['stage']
            kwargs.pop('stage')
        else: stage= None
        if kwargs.has_key('d1'):
            d1= kwargs['d1']
            kwargs.pop('d1')
        else: d1= self._filters[0]+'-'+self._filters[1]
        if kwargs.has_key('d2'):
            d2= kwargs['d2']
            kwargs.pop('d2')
        else: d2= self._filters[0]
        if kwargs.has_key('ignore_gaps'):
            ignore_gaps= kwargs['ignore_gaps']
            kwargs.pop('ignore_gaps')
        else:
            ignore_gaps= False
        #get isochrone
        try:
            iso= self(logage,Z=Z,feh=feh,afe=afe,maxm=maxm,stage=stage)
        except IOError:
            if ignore_gaps: return None
            else: 
                raise IOError("No isochrone found for this logage/metallicity combination\nUse ignore_gaps=True to ignore non-existant isochrones")
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
                kwargs['yrange']= [numpy.amax(y)+0.3,numpy.amin(y)-0.3]
        #plot
        return bovy_plot.bovy_plot(x,y,*args,**kwargs)

def Z2FEH(z,zsolar=None,parsec=False):
    """Convert Z to FeH assuming zsolar"""
    if parsec:
        if zsolar is None: zsolar= 0.0152
        return numpy.log10(z/(1.-0.2485-2.78*z))-math.log10(zsolar/(1.-0.2485-2.78*zsolar))
    else:
        if zsolar is None: zsolar= _ZSOLAR
        return numpy.log10(z)-math.log10(zsolar)

def FEH2Z(feh,zsolar=None,parsec=False):
    """Convert FeH to Z assuming zsolar"""
    if parsec:
        if zsolar is None: zsolar= 0.0152
        zx= 10.**(feh+math.log10(zsolar/(1.-0.2485-2.78*zsolar)))
        return (zx-0.2485*zx)/(2.78*zx+1.)
    else:
        if zsolar is None: zsolar= _ZSOLAR
        return 10.**(feh+math.log10(zsolar))

def logg(logL,logTe,mass):
    """
    NAME:
       logg
    PURPOSE:
       calculate log g from luminosity, teff, and mass
    INPUT:
       logL - log Luminosity (solar units)
       logTe- log effective temperature (log K)
       mass - mass (solar units)
    OUTPUT:
       log g (log cm/s^2)
    HISTORY:
       2012-08-16 - Written - Bovy (IAS)
    """
    logR= -2.*(logTe-_LOGTESUN)+0.5*logL
    return numpy.log10(mass)-2.*logR+_LOGGSUN

def dict2recarray(dict):
    nEntries= len(dict.keys())
    nOut= len(dict[dict.keys()[0]])
    out= numpy.zeros(nOut,dtype={'names':dict.keys(),
                                 'formats':[numpy.float64 for ii in range(nEntries)]})
    for ii in range(nEntries):
        out[dict.keys()[ii]]= dict[dict.keys()[ii]]
    return out.view(numpy.recarray)
