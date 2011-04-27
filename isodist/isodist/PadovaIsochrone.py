import csv
import numpy as nu
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
       read_padova_isochrone('../data/isoc_z0.dat',filters=['J','H','Ks'])
    HISTORY:
       2011-04-26 - Written - Bovy (NYU)
    """
    dialect= csv.excel
    dialect.skipinitialspace=True
    reader= csv.reader(open(name,'r'),delimiter='\t',
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
    Flum= []
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
        Flum.append(float(row[8+nfilters]))
    #Load everything into a dictionary
    outDict= {}
    outDict['logage']= nu.array(logage)
    outDict['M_ini']= nu.array(M_ini)
    outDict['M_act']= nu.array(M_act)
    outDict['logL']= nu.array(logL)
    outDict['logTe']= nu.array('logTe')
    outDict['logg']= nu.array(logg)
    outDict['mbol']= nu.array(mbol)
    outDict['Flum']= nu.array(Flum)
    for ii in range(nfilters):
        thismag= []
        for jj in range(len(mags)):
            thismag.append(mags[jj][ii])
        outDict[filters[ii]]= nu.array(thismag)
    return outDict
