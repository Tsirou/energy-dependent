import sys

import names
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord


#==============================================================================
# Functions
#==============================================================================
    
# In order to have an output countdown
def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s --- %s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()


# Gnomonic projection for computing the distance between two objects
def angular_separation(l1, b1, l2, b2):
    aa                    = np.sin(b1) * np.sin(b2) + np.cos(b1) * np.cos(b2) * np.cos(l2 - l1)
    angle                 = np.arccos(aa)
    return angle
#==============================================================================

# Conversion from pixels to celestial coordinates using the center position and the resolution of the pic
def convert_image(cx,cy,cim,xpos,ypos,res):
    cel  = cx - (xpos - cim) * res / (np.cos(cy * np.pi / 180.))
    esti = cy + (ypos - cim) * res

    return cel,esti

# Conversion from pixels to world celestial coordinates using the FITS header and astropy python library
def convert_pix2wcs(fits_filename,x,y):
    header         = fits.getheader(fits_filename)
    w              = wcs.WCS(header)
    
    cel1, cel2   = w.wcs_pix2world(x, y,names.random_convert_shift)

    # if(names.analysis == 0):
    #     c = SkyCoord(l=cel1*u.degree, b=cel2*u.degree, frame='galactic')
    #
    #     print "\nalpha = ",c.fk5.ra,"\ndelta = ",c.fk5.dec,"deg\n LON,LAT = ",cel1,cel2,"\n"
    #
    #     ra  = c.fk5.ra.deg
    #     dec = c.fk5.dec.deg
    # else:
    ra  = cel1
    dec = cel2

    return ra,dec

# For the error conversion from pixels to world celestial coordinates    
def error_convert_ra(fits_filename,x):
    header  = fits.getheader(fits_filename)
    w       = wcs.WCS(header)

    ref     = w.wcs_pix2world(0,0,names.random_convert_shift)
    deg     = w.wcs_pix2world(x,0,names.random_convert_shift)
    # returns value in seconds
    return abs(deg[0]-ref[0])*12*3600/180.      

def error_convert_dec(fits_filename,x):
    header  = fits.getheader(fits_filename)
    w       = wcs.WCS(header)

    ref     = w.wcs_pix2world(0,0,names.random_convert_shift)
    deg     = w.wcs_pix2world(0,x,names.random_convert_shift)
    #returns value in arcsec
    return abs(deg[1]-ref[1])*3600
 
#==============================================================================    
# Converting RA,DEC from degrees to RA in h:m:s and DEC deg:arcmin:arcsec [display]
def convert_deg2timearcs(ra,dec):
    
    ra_total_hour     = ra * 24. / 360.
    ra_hour           = int(ra_total_hour)
    ra_min            = int((ra_total_hour - ra_hour) * 60. )
    ra_sec            = np.round((ra_total_hour * 3600. - (ra_hour * 3600.) - (ra_min * 60.)), 1)
    
    dec_deg           = int(dec)
    dec_total_arcsec  = (abs(dec-dec_deg))* np.pi / 180. / 4.848e-6
    dec_arcmin        = int(dec_total_arcsec) / 60
    dec_arcsec        = np.round(dec_total_arcsec - dec_arcmin*60.,1)
    
    print "\nRA  = ",ra_hour,"h ",ra_min,"min ",ra_sec,"sec\n"
    print "\nDEC = ",dec_deg,"deg ", dec_arcmin,"arcmin ", dec_arcsec," arcsec\n"
     
#==============================================================================
# Rebinning functions  
def rebin( a, newshape ):
    assert len(a.shape) == len(newshape)
    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]

             
def rebin_factor( a, newshape ):
        #newshape must be a factor of a.shape.
        assert len(a.shape) == len(newshape)
        assert not np.sometrue(np.mod( a.shape, newshape ))

        slices = [ slice(None,None, old/new) for old,new in zip(a.shape,newshape) ]
        return a[slices]

def rebin_med(a,factor):

    b_i = len(a)/factor
    b_j = len(a)/factor
    
    b   = np.zeros((b_i,b_j))
    
    for ii in range(0,b_i):
        for jj in range(0,b_j):
            
            for i in range(ii*factor,(ii+1)*factor):
                for j in range(jj*factor,(jj+1)*factor):
                    
                    b[ii][jj]  = b[ii][jj] + a[i][j]

            b[ii][jj] = b[ii][jj] / (factor**2)  
    return b


def rebin_med_non_avg(a, factor):
    b_i = len(a) / factor
    b_j = len(a) / factor

    b = np.zeros((b_i, b_j))

    for ii in range(0, b_i):
        for jj in range(0, b_j):

            for i in range(ii * factor, (ii + 1) * factor):
                for j in range(jj * factor, (jj + 1) * factor):
                    b[ii][jj] = b[ii][jj] + a[i][j]

            b[ii][jj] = b[ii][jj]
    return b

    
def rebin_red(a,factor):

    b_i = len(a) * factor
    b_j = len(a) * factor
    
    b   = np.zeros((b_i,b_j))
    

    for i in range(0, len(a)):
        for j in range(0, len(a)):

            for ii in range(i * factor, (i+1)* factor):
                for jj in range(j * factor,(j+1) * factor):
                           
                    b[ii][jj]  = b[ii][jj] + a[i][j]

            b[ii][jj] = b[ii][jj] #/ (factor**2)  
    return b    
    
    
#==============================================================================
# Cropping function    
def crop(a,new_dim):
    
    dim_res = (len(a) - new_dim) / 2
    b       = a[dim_res:len(a)-dim_res,dim_res:len(a)-dim_res]
    return b


# ==============================================================================
# Cicular cropping function
def crop_circle(a, rad):
    circle    = np.zeros((len(a),len(a)))
    for i in range(0,len(a)):
        for j in range(0,len(a)):
            if( (i)**2 + (j)**2 <= rad**2):
                circle[len(a)/2 - 1 + i, len(a)/2 - 1 + j] = 1
                circle[len(a)/2 - 1 - i, len(a)/2 - 1 + j] = 1
                circle[len(a)/2 - 1 + i, len(a)/2 - 1 - j] = 1
                circle[len(a)/2 - 1 - i, len(a)/2 - 1 - j] = 1

    for i in range(0,len(a)):
        for j in range(0,len(a)):
            if(circle[i,j] == 0):
                a[i,j] = 0
    b         = a[:,:]

    return b

#==============================================================================
# Filling an image to an extent with zeros    
def fill_fits(a,dim_fill):
    b = np.zeros((dim_fill,dim_fill))
    if (dim_fill < len(a)):
        print "Cannot fill!\n"
    for i in range(dim_fill/2-len(a)/2 ,dim_fill/2+len(a)/2):
        for j in range(dim_fill/2-len(a)/2 ,dim_fill/2+len(a)/2):
            b[i][j] = a[i-(dim_fill/2-len(a)/2)][j-(dim_fill/2-len(a)/2)]
    return b

#==============================================================================
# Normalizing a fit image with a given norm
def norm_fits(a,norm):
    
    dim_im = len(a)
    b      = np.zeros((dim_im,dim_im))
    
    for i in range(0,dim_im):
        for j in range(0,dim_im):
           b[i,j] = a[i,j] * norm 

    return b
            
#==============================================================================
# Converting the FWHM from pixels to sigma in arcmin
def fwhm_sigma( fwhm, factor):
    size  = fwhm / (2. * np.sqrt(2. * np.log(2.)))    
    size  = size * factor * 60.
    return size
    
#==============================================================================
# Finding a given percentage of the emission inside a radius
def radius_percent(percent,hist,pix2degree):
    all_sum = 0
    for i in range(0,len(hist)):
        progress(i,len(hist), '')
        for j in range(0,len(hist)):
            all_sum  = all_sum + hist[i][j]
    
    sum_68 = hist[int(len(hist)/2.)][int(len(hist)/2.)]
   
    for k in range(1,int(len(hist)/2.)):
        if(sum_68 < all_sum * percent / 100.): 
            for l in range(-k,k+1):
                sum_68 = sum_68 + \
hist[int(len(hist)/2.) + l][int(len(hist)/2.) + k] + \
hist[int(len(hist)/2.) + l][int(len(hist)/2.) - k] + \
hist[int(len(hist)/2.) - k][int(len(hist)/2.) + l] + \
hist[int(len(hist)/2.) + k][int(len(hist)/2.) + l]  
            sum_68 = sum_68 - \
hist[int(len(hist)/2.) + k][int(len(hist)/2.) + k] - \
hist[int(len(hist)/2.) - k][int(len(hist)/2.) + k] - \
hist[int(len(hist)/2.) + k][int(len(hist)/2.) - k] - \
hist[int(len(hist)/2.) - k][int(len(hist)/2.) - k] 

            radius_it  = k * pix2degree * 60.
     
    return all_sum,sum_68,radius_it 
            
    
#==============================================================================
# Computing the mean and the standart deviation
def stat(x):    
    x2   = []
    
    for element in x:
       x2.append(element**2.)

    sigma  = np.sqrt(np.mean(x2) - (np.mean(x))**2.)
    mu     =  np.mean(x)
        
    return mu,sigma

def log_likelihood_gauss(x,mu,sigma):
    l  = (-0.5 * len(x) * np.log(2.*np.pi*sigma**2)) - ((1./(2.*sigma**2)) * sum((x-mu)**2))
    return l    

def AIC(l,k):
    return -2. * l + 2.*k    

#==============================================================================
# Definition of gaussians
def gaussian(x, x0, norm, sigma):
    a  = norm * 1./(sigma * np.sqrt(2. * np.pi))
    return a * np.exp(-((x - x0)**2 / (2 * sigma**2)))

def lorentz(x, x0, norm, gamma):
    a  = norm * 2. / (np.pi * gamma)
    return a * (1./ ( 1 + ((x - x0)/ (gamma / 2))**2))             

def lorentz_gaussian(x,x0,norm_lor,norm_gauss,sigma,gamma): 
    a  = norm_lor * 2. / (np.pi * gamma)
    b  = norm_gauss * 1./(sigma * np.sqrt(2. * np.pi))
    return  (a * (1./ ( 1 + ((x - x0)/ (gamma / 2))**2)) )+(b * np.exp(-((x - x0)**2 / (2 * sigma**2))))          

def gaussian_gaussian(x,x0,norm_gauss1,norm_gauss2,sigma1,sigma2): 
    a  = norm_gauss1 * 1./(sigma1 * np.sqrt(2. * np.pi))
    b  = norm_gauss2 * 1./(sigma2 * np.sqrt(2. * np.pi))
    return  (a * np.exp(-((x - x0)**2 / (2 * sigma1**2))) )+(b * np.exp(-((x - x0)**2 / (2 * sigma2**2))))  

def triple_gaussian(x,x0,norm_gauss1,norm_gauss2,norm_gauss3,sigma1,sigma2,sigma3): 
    a      = norm_gauss1 * 1./(sigma1 * np.sqrt(2. * np.pi))
    b      = norm_gauss2 * 1./(sigma2 * np.sqrt(2. * np.pi))
    c      = norm_gauss3 * 1./(sigma3 * np.sqrt(2. * np.pi))

    triple = (a * np.exp(-((x - x0)**2 / (2 * sigma1**2))) ) + \
(b * np.exp(-((x - x0)**2 / (2 * sigma2**2)))) + \
(c * np.exp(-((x - x0)**2 / (2 * sigma3**2))))
    return triple


def affine(x,x1,y1,x2,y2):

    a   = (y1 - y2) / (x1 - x2)
    b   = y1 - a * x1
    if( abs( b - (y2 - a * x2)) > 1.0e-2):
        print "\n Problem with this linear function, check input coordinates!\n"

    return a * x + b

# Found on the web
def radial_profile(data, center):
    x, y = np.indices((data.shape))
    r    = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r    = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr   = np.bincount(r.ravel())
    rp   = tbin / nr
    return rp

def factorial(n):
    x = 1
    for i in range(1, (n+1)):
        x = x * i
    return x

def Hermite_polynoms(x,n):
    sum_H   = 0

    if( n % 2 == 0):
        for l in range(0,(n/2) + 1):
            sum_H = sum_H + ( (-1.)**( (n/2) - l)  / (factorial(2 * l) * factorial((n/2) - l)) * (2 * x)**(2 * l))
        H = factorial(n) * sum_H
    if( n % 2 == 1):
        for l in range(0,((n - 1)/2) + 1):
            sum_H = sum_H + ((-1.)**( ((n - 1)/2) - l) / (factorial(2 * l + 1) * factorial( int( (n-1)/2) - l)) * (2 * x)**(2 * l + 1))
        H = factorial(n) * sum_H

    return H


def pixarea_correction(exposure_fits, ra, dec, resolution, corr_exposure_fits):

    dimension            = len(exposure_fits)

    for i in range(0, dimension):
        corr_exposure_fits[:][i]   = np.cos( (-1.0 * dec - (i - (dimension/2. - 0.5)) * resolution) / 180. * np.pi )

    corr_exposure_fits[:][:] = corr_exposure_fits[:][:] / (np.cos( (-1.0 * dec) / 180. * np.pi ))

    corr_exposure_fits[:][:] = corr_exposure_fits[:][:] * exposure_fits[:][:]

    return corr_exposure_fits