from mechanize import Browser
from math import exp, log
import os
import scipy as sc
import time
isodir= 'isochrones/'
def get_isochrones(nzs=301,ages=[6.6,10.2,0.012],wait=0.):
    zs= sc.linspace(0.0001,0.03,nzs)
    br= Browser()
    for z in zs:
        zstring= str(z)
        br.open('http://stev.oapd.inaf.it/cgi-bin/cmd')
        for form in br.forms():#There is only one, hopefully!
            br.form= form
            br["photsys_file"]=["tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat"]
            #br["photsys_file"]=["tab_mag_odfnew/tab_mag_ubvrijhk.dat"]
            br["isoc_val"]= ["1"]
            br["isoc_zeta0"]= zstring
            br["isoc_lage0"]=str(ages[0])
            br["isoc_lage1"]=str(ages[1])
            br["isoc_dlage"]=str(ages[2])
            br.submit()
        link=br.find_link()
        filename= link.text
        savefilename='Z-'+zstring+'.dat.gz'
        os.system('wget http://stev.oapd.inaf.it/~lgirardi/tmp/'+filename)
        os.system('mv '+filename+' '+savefilename)
        os.system('mv '+savefilename+' '+isodir)
        if wait != 0:
            time.sleep(wait)
    return 0

if __name__ == '__main__':
    get_isochrones(301,wait=5.)
