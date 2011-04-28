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
