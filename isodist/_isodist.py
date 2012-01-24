import numpy as nu
from scipy.maxentropy import logsumexp
from Isochrone import Isochrone
from PadovaIsochrone import PadovaIsochrone
def eval_distpdf(ds,mdict=None,mivardict=None,logg=None,logg_ivar=None,
                 teff=None,teff_ivar=None,logage=None,logage_ivar=None,
                 Z=None,Z_ivar=None,feh=None,feh_ivar=None,
                 afe=None,afe_ivar=None,
                 padova=None,padova_type=None,
                 normalize=False,
                 ageprior=None):
    """
    NAME:
       eval_distpdf
    PURPOSE:
       evaluate the distance PDF for an object
    INPUT:
       ds- list or ndarray of distance (or a single distance), in kpc
       mdict= dictionary of apparent magnitudes (e.g., {'J':12.,'Ks':13.})
       mivardict= dictionary of magnitude inverse variances (matched to mdict)
       logg= observed logg
       logg_ivar= inverse variance of logg measurement
       teff= observed T_eff [K]
       logg_ivar= inverse variance of T_eff measurement
       logage= observed log_10 age [Gyr]
       logage_ivar= inverse variance of log_10 age measurement
       Z= observed metallicity
       Z_ivar= inverse variance of Z measurement
       feh= observed metallicity (alternative to Z)
       feh_ivar= inverse variance of FeH measurement
       afe= observed [\alpha/Fe]
       afe_ivar= [\alpha/Fe] inverse variance
       padova= if True, use Padova isochrones, 
               if set to a PadovaIsochrone objects, use this
       padova_type= type of PadovaIsochrone to use (e.g., 2mass-spitzer-wise)
       normalize= if True, normalize output PDF (default: False)
       ageprior= - None: flat in log age
                 - flat: flat in age
    OUTPUT:
       log of probability
    HISTORY:
       2011-04-28 - Written - Bovy (NYU)
    """
    #load isochrones
    if not padova is None and isinstance(padova,PadovaIsochrone):
        iso= padova
    elif not padova is None and isinstance(padova,bool) and padova:
        iso= PadovaIsochrone(type=padova_type)
    #Parse metallicity info
    if not feh is None: raise NotImplementedError("'feh' not yet implemented")
    #set up output
    if isinstance(ds,(list,nu.ndarray)):
        scalarOut= False
        if isinstance(ds,list):
            _ds= nu.array(ds)
        else: _ds= ds
    elif isinstance(ds,float):
        scalarOut= True
        _ds= [ds]
    #Pre-calculate all absolute magnitudes
    absmagdict= {}
    for key in mdict.keys():
        absmagdict[key]= -_distmodulus(_ds)+mdict[key]        
    #loop through isochrones
    ZS= iso.Zs()
    logages= iso.logages()
    allout= nu.zeros((len(_ds),len(ZS),len(logages)))
    for zz in range(len(ZS)):
        for aa in range(len(logages)):
            thisiso= iso(logages[aa],Z=ZS[zz])
            dmpm= nu.roll(thisiso['M_ini'],-1)-thisiso['M_ini']
            loglike= nu.zeros((len(_ds),len(thisiso['M_ini'])-1))
            loglike-= nu.log(thisiso['M_ini'][-1])
            for ii in range(1,len(thisiso['M_ini'])-1):
                if dmpm[ii] > 0.: 
                    loglike[:,ii]+= nu.log(dmpm[ii])
                else: 
                    loglike[:,ii]= nu.finfo(nu.dtype(nu.float64)).min
                    continue #no use in continuing here
                if not teff is None:
                    loglike[:,ii]-= (teff-10**thisiso['logTe'][ii])**2.*teff_ivar
                if not logg is None:
                    loglike[:,ii]-= (logg-thisiso['logg'][ii])**2.*logg_ivar
                for key in mdict.keys():
                    #print absmagdict[key][2], thisiso[key][ii]
                    loglike[:,ii]-= (absmagdict[key]-thisiso[key][ii])**2.\
                        *mivardict[key]
            #marginalize over mass
            for jj in range(len(_ds)):
                allout[jj,zz,aa]= logsumexp(loglike[jj,:])
            #add age constraint and prior
            if not logage is None:
                allout[:,zz,aa]+= -(logage-logages[aa])**2.*logage_ivar
            if not ageprior is None:
                if isinstance(ageprior,str) and ageprior.lower() == 'flat':
                    allout[:,zz,aa]+= logages[aa]
        #add Z constraint and prior
        if not Z is None:
            allout[:,zz,:]+= -(Z-ZS[zz])**2.*Z_ivar
    #prepare final output
    out= nu.zeros(len(_ds))
    for jj in range(len(_ds)):
        out[jj]= logsumexp(allout[jj,:,:])
    if normalize:
        out-= logsumexp(out)
    #return
    if scalarOut: return out[0]
    else: return out

def _distmodulus(d):
    return 5.*nu.log10(d/.01)
