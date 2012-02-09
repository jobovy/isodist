###############################################################################
# isodist.imf: various IMF models
#
#  These all return dN/dM
###############################################################################
import numpy
from scipy import integrate
_LOGMOLOGNORMALCHABRIER2001= numpy.log10(0.1)
_S2LOGNORMALCHABRIER2001= 0.627**2.
_LOGLN= numpy.log(10.)
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

def _intLognormalChabrier2001Integrand(m):
    return lognormalChabrier2001(m,int=False)
