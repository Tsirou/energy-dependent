import sys
import conversion
import numpy as np
from astropy.coordinates import SkyCoord


analysis   = ''

ra_psr     = 228.5292
dec_psr    = -59.1575
glon_psr   = 320.34
glat_psr   = -1.20

xray_ampl  = 82990.7057651
comp_ampl  = 184.597258064
comp_xpos  = 198.162464344
comp_ypos  = 197.955058846
comp_ellip = 0.
comp_theta = 0.
comp_fwhm  = 27.0503797309




sys.stdout.write("hgps or pa?\n")
analysis = sys.stdin.readline()




if(analysis.find("hgps") != -1):
    binsize    = 0.02
    pic_center = 100.5

    coord1, coord2 = conversion.convert_image(glon_psr,glat_psr,pic_center,comp_xpos,comp_ypos,binsize)

    equ = SkyCoord(l=coord1,b=coord2,frame='galactic',unit='deg').fk5

    print "Component amplitude / X-ray amplitude : ",comp_ampl/xray_ampl

    print "l : ", coord1, " , b : ", coord2

    conversion.convert_deg2timearcs(equ.ra.deg,equ.dec.deg)

    print "Major axis :",conversion.fwhm_sigma(comp_fwhm,binsize),"\n"

    print "Minor axis :",conversion.fwhm_sigma(comp_fwhm*(1-comp_ellip),binsize),"\n"

    if (comp_theta >= 0.):
        print "Inclination : ",comp_theta * 180. / np.pi,"\n"
    if (comp_theta < 0.):
        print "Inclination : ",( 360. + (comp_theta * 180. / np.pi)),"\n"


if(analysis.find("pa") != -1):
    binsize    = 0.01
    pic_center = 200.5

    coord1, coord2 = conversion.convert_image(ra_psr,dec_psr,pic_center,comp_xpos,comp_ypos,binsize)

    print coord1,coord2

    equ = SkyCoord(ra=coord1,dec=coord2,unit='deg')

    print "Component amplitude / X-ray amplitude : ", comp_ampl / xray_ampl

    conversion.convert_deg2timearcs(equ.ra.deg,equ.dec.deg)

    print "l : ", equ.galactic.l.deg, " , b : ", equ.galactic.b.deg

    print "Major axis :",conversion.fwhm_sigma(comp_fwhm,binsize),"\n"

    print "Minor axis :",conversion.fwhm_sigma(comp_fwhm*(1-comp_ellip),binsize),"\n"

    if (comp_theta >= 0.):
        print "Inclination : ",comp_theta * 180. / np.pi,"\n"
    if (comp_theta < 0.):
        print "Inclination : ",( 360. + (comp_theta * 180. / np.pi)),"\n"





