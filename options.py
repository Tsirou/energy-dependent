import sys
import names

#Save option
def save_fits():
    sys.stdout.write("\nDo you want to save the results? y/n \n")
    save_option = names.save_option
    while(save_option.find('\n') == -1):
        save_option  = sys.stdin.readline()
    if (save_option != "y\n" and save_option != "n\n"):
        sys.stdout.write("\n Please choose y or n ...\n")

    return save_option

    
# Selecting a kind of fit for the gamma-like candidate image (either a simple gaussian fit or mine)
def simple():
    sys.stdout.write("\nWhat kind of fit? simple/complex \n")
    fit_option = names.fit_option
    while(fit_option.find('\n') == -1):
        fit_option  = sys.stdin.readline()
    if (fit_option != "simple\n" and fit_option != "complex\n"):
        sys.stdout.write("\n Please choose simple or complex keywords.\n")
    sys.stdout.write("\n\n")
    
    return fit_option
 
#Loading an HESS PSF image    
def use_psf():   
    sys.stdout.write("\nUse of the HESS PSF? y/n \n")
    psf                            = sys.stdin.readline()

    return psf
    
def nb_gauss():
    sys.stdout.write("\nFor fitting the source, how many gaussians? (G1, G2 or G3) :\n")
    psf_gauss  = sys.stdin.readline()     

    return psf_gauss
    
def psf_fit():
    sys.stdout.write("\nWould you like to fit the given PSF? y/n :\n")
    fit_psf  = sys.stdin.readline()     

    return fit_psf    

def cash_profile():
    sys.stdout.write("\nCompute something similar to a likelihood profile? y/n :\n")
    l_p  = sys.stdin.readline()     

    return l_p 

def bkg_model():
    sys.stdout.write("\nWhich background model? temp (tmp)/gamma-like? (gam) or (rng) :\n")
    bkg_md  = sys.stdin.readline()     

    return bkg_md         