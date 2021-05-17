###############################################################################
# isodist.imf: various IMF models
#
#  These all return dN/dM
###############################################################################
import numpy
from scipy import integrate
_LOGMOLOGNORMALCHABRIER2001= numpy.log10(0.1)
_S2LOGNORMALCHABRIER2001= 0.627**2.
_LOGMOLOGNORMALCHABRIER2003= numpy.log10(0.079)
_S2LOGNORMALCHABRIER2003= 0.69**2.
_ALPHAEXPONENTIALCHABRIER2001= 3.3
_BETAEXPONENTIALCHABRIER2001= 0.25
_MOEXPONENTIALCHABRIER2001= 716.4
_LOGLN= numpy.log(10.)
import sys
_PY3= sys.version > '3'
if _PY3:
    long= int
def lognormalChabrier2001(m,int=False):
    """
    NAME:
       lognormalChabrier2001
    PURPOSE:
       IMF2 from Chabrier (2001), ApJ, 554, 1274
    INPUT:
       m - mass in solar masses
       int= (default: False) if True, return integrated N(<m)
    OUTPUT:
       dN/dm
    HISTORY:
       2012-02-08 - Written - Bovy (IAS)
    """
    if int:
        if isinstance(m,(long,float)):
            m= [m]
            scalarOut= True
        else:
            scalarOut= False
        out= numpy.zeros(len(m))
        for ii in range(len(m)):
            out[ii]= integrate.quad(_intLognormalChabrier2001Integrand,
                                    0.,
                                    m[ii])[0]
        if isinstance(m,list): return list(out)
        elif scalarOut: return out[0]
        else: return out
    else:
        return 0.141/m*numpy.exp(-(numpy.log10(m)-_LOGMOLOGNORMALCHABRIER2001)**2./2./_S2LOGNORMALCHABRIER2001)/_LOGLN

def exponentialChabrier2001(m,int=False):
    """
    NAME:
       exponentialChabrier2001
    PURPOSE:
       IMF3 from Chabrier (2001), ApJ, 554, 1274
    INPUT:
       m - mass in solar masses
       int= (default: False) if True, return integrated N(<m)
    OUTPUT:
       dN/dm
    HISTORY:
       2012-02-08 - Written - Bovy (IAS)
    """
    if int:
        if isinstance(m,(long,float)):
            m= [m]
            scalarOut= True
        else:
            scalarOut= False
        out= numpy.zeros(len(m))
        for ii in range(len(m)):
            out[ii]= integrate.quad(_intExponentialChabrier2001Integrand,
                                    0.,
                                    m[ii])[0]
        if isinstance(m,list): return list(out)
        elif scalarOut: return out[0]
        else: return out
    else:
        return 3.*m**-_ALPHAEXPONENTIALCHABRIER2001*numpy.exp(-(_MOEXPONENTIALCHABRIER2001/m)**_BETAEXPONENTIALCHABRIER2001)

def kroupa2003(m,int=False):
    """
    NAME:
       kroupa2003
    PURPOSE:
       ''universal'' IMF from Kroupa (2003), MNRAS 322, 231
    INPUT:
       m - mass in solar masses
       int= (default: False) if True, return integrated N(<m)
    OUTPUT:
       dN/dm
    HISTORY:
       2012-02-08 - Written - Bovy (IAS)
    """
    if int:
        if isinstance(m,(long,float)):
            m= [m]
            scalarOut= True
        else:
            scalarOut= False
        out= numpy.zeros(len(m))
        for ii in range(len(m)):
            out[ii]= integrate.quad(_intKroupa2003Integrand,
                                    0.,
                                    m[ii])[0]
        if isinstance(m,list): return list(out)
        elif scalarOut: return out[0]
        else: return out
    else:
        if isinstance(m,(long,float)):
            if m < 0.08: return (m/0.08)**-0.3
            elif m < 0.5: return (m/0.08)**-1.3
            else: return (m/0.5)**-2.3*(0.5/0.08)**-1.3
        elif isinstance(m,list):
            m= numpy.array(m)
            return kroupa2003(m)
        elif isinstance(m,numpy.ndarray):
            out= numpy.zeros(len(m))
            out[(m < 0.08)]= (m[(m < 0.08)]/0.08)**-0.3
            out[(m >= 0.08)*(m < 0.5)]= (m[(m >= 0.08)*(m < 0.5)]/0.08)**-1.3
            out[(m >= 0.5)]= (m[(m >= 0.5)]/0.5)**-2.3*(0.5/0.08)**-1.3
            return out

def chabrier2003(m,int=False):
    """
    NAME:
       chabrier2003
    PURPOSE:
       IMF from Table 1 from Chabrier (2003), PASP, 115, 763
    INPUT:
       m - mass in solar masses
       int= (default: False) if True, return integrated N(<m)
    OUTPUT:
       dN/dm
    HISTORY:
       2012-02-08 - Written - Bovy (IAS)
    """
    if int:
        if isinstance(m,(long,float)):
            m= [m]
            scalarOut= True
        else:
            scalarOut= False
        out= numpy.zeros(len(m))
        for ii in range(len(m)):
            out[ii]= integrate.quad(_intChabrier2003Integrand,
                                    0.,
                                    m[ii])[0]
        if isinstance(m,list): return list(out)
        elif scalarOut: return out[0]
        else: return out
    else:
        if isinstance(m,(long,float)):
            if m < 1.: return 0.158/m*numpy.exp(-(numpy.log10(m)-_LOGMOLOGNORMALCHABRIER2003)**2./2./_S2LOGNORMALCHABRIER2003)/_LOGLN
            else: return m**-2.3*0.158*numpy.exp(-_LOGMOLOGNORMALCHABRIER2003**2./2./_S2LOGNORMALCHABRIER2003)/_LOGLN
        elif isinstance(m,list):
            m= numpy.array(m)
            return chabrier2003(m)
        elif isinstance(m,numpy.ndarray):
            out= numpy.zeros(len(m))
            out[(m < 1.)]= 0.158/m[(m < 1.)]*numpy.exp(-(numpy.log10(m[(m < 1.)])-_LOGMOLOGNORMALCHABRIER2003)**2./2./_S2LOGNORMALCHABRIER2003)/_LOGLN
            out[(m >= 1.)]= m[(m >= 1.)]**-2.3*0.158*numpy.exp(-_LOGMOLOGNORMALCHABRIER2003**2./2./_S2LOGNORMALCHABRIER2003)/_LOGLN
            return out

def _intLognormalChabrier2001Integrand(m):
    return lognormalChabrier2001(m,int=False)
def _intExponentialChabrier2001Integrand(m):
    return exponentialChabrier2001(m,int=False)
def _intKroupa2003Integrand(m):
    return kroupa2003(m,int=False)
def _intChabrier2003Integrand(m):
    return chabrier2003(m,int=False)
