from Isochrone import Isochrone
def eval_distpdf(ds,mdict=None,mivardict=None,logg=None,logg_ivar=None,
                 teff=None,teff_ivar=None,age=None,age_ivar=None,
                 Z=None,Z_ivar=None,feh=None,feh_ivar=None,
                 afe=None,afe_ivar=None,
                 padova=None,padova_type=None):
    """
    NAME:
       eval_distpdf
    PURPOSE:
       evaluate the distance PDF for an object
    INPUT:
       ds- list or ndarray of distance (or a single distance)
       mdict= dictionary of apparent magnitudes (e.g., {'J':12.,'Ks':13.})
       mivardict= dictionary of magnitude inverse variances (matched to mdict)
       logg= observed logg
       logg_ivar= inverse variance of logg measurement
       teff= observed T_eff
       logg_ivar= inverse variance of T_eff measurement
       age= observed age
       age_ivar= inverse variance of age measurement
       Z= observed metallicity
       Z_ivar= inverse variance of Z measurement
       feh= observed metallicity (alternative to Z)
       feh_ivar= inverse variance of FeH measurement
       afe= observed [\alpha/Fe]
       afe_ivar= [\alpha/Fe] inverse variance
       padova= if True, use Padova isochrones, 
               if set to a PadovaIsochrone objects, use this
       padova_type= type of PadovaIsochrone to use (e.g., 2mass-spitzer-wise)
    OUTPUT:
       log of probability
    HISTORY:
       2011-04-28 - Written - Bovy (NYU)
    """
    return None
