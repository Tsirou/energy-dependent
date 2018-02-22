import sys,os,string
import conversion
import names 

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord


def display_results(save_fit,psf,psf_gauss,bkg_option):
    
    if (save_fit.find("yes") != -1):
        nb_displays = 2
    else :
        nb_displays = 1

    if (save_fit.find("msh") != -1):
        saved_results = names.save_path[names.analysis] + "conf_"+psf_gauss[:2]+".txt"
        saved_fitinfo = names.save_path[names.analysis] + "fit_info_"+psf_gauss[:2]+".out"
    if (save_fit.find("ls") != -1):
        saved_results = names.save_path + "conf_ls5039"+psf_gauss[:2]+".txt"
        # Not implemented saved_fitinfo = names.save_path + "fit_info_ls5039_"+psf_gauss[:2]+".out"
    if (save_fit.find("X") != -1):
        saved_results = names.save_path[names.analysis] + "conf_"+psf_gauss[:3]+".txt"
        saved_fitinfo = names.save_path[names.analysis] + "fit_info_"+psf_gauss[:3]+".out"


    for disp in range(0,nb_displays):
        if(nb_displays == 1):
            temporary     = sys.stdout
        if(nb_displays == 2):
            temporary     = sys.stdout
            sys.stdout = open(names.results_path + "final_fit_results"+psf_gauss[:2]+"_"+psf[:1]+"_"+bkg_option[:3]+".txt", 'w')
        

    if(psf.find("n") != -1):
        fits_filename              = names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5]+"_"+psf_gauss[:2]+bkg_option[:2]+"_model.fits"
    if(psf.find("y") != -1):
        if( save_fit.find("msh") != -1 or save_fit.find("ls") != -1):
            fits_filename              = names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:2] + "_model.fits"
        if( save_fit.find("X") != -1):
            print("Not implemented ...")
            fits_filename = names.path_template[names.analysis] + names.filename_gamma[names.analysis]

    hdu                            = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])
    factor_pix2deg_gamma           = abs(hdu[names.hdu_ext].header['cdelt1'])
    hdu.close()

    print("\n ***********\n **********\n ________________________\n ####### Results ######\n")
    print("Object      : " + names.source[names.analysis] + "\n")
    print("Type of fit :\n   PSF model :"+psf[0]+"\n   Gaussians : "+psf_gauss[1:2]+"\n   Bkg model : "+bkg_option)

    nf,cash_stats    =  np.loadtxt(saved_fitinfo,dtype=float,usecols=(0,1),unpack=True,skiprows=1)
    cash             = cash_stats[-1] 

    with file(saved_fitinfo) as f:
        one_row = f.readline()
    dof = len(one_row.split()) - 3
    f.close()

    if (save_fit.find("X") != -1):
        print "\n STATISTICS :\n   CASH  = ", cash, "\n   d.o.f : ", dof, "\n   AIC   = ", cash + 2. * (dof + 1)
    else:
        print "\n STATISTICS :\n   CASH  = ",cash,"\n   d.o.f : ",dof,"\n   AIC   = ",cash + 2.*dof

    conf_names,equals    =  np.loadtxt(saved_results,dtype=str,usecols=(0,1),unpack=True,skiprows=1)
    column               = ["parnames","parvals","parmins","parmaxes"]

    if (save_fit.find("X") != -1):
        if(psf_gauss.find('G') == -1 and psf_gauss.find('sh') == -1 and psf_gauss.find('dc') == -1):
            param_names          = ["ampl"]

        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            param_names          = ["fwhm","xpos","ypos","ellip","theta","xrays_mod.ampl","gcomp.ampl"]
        if(psf_gauss.find("2G") != -1):
            param_names          = ["fwhm","xpos","ypos","xrays_mod.ampl","gcomp.ampl"]
        if(psf_gauss.find("sh") != -1):
            param_names          = ["r0","xpos","ypos","width","xrays_mod.ampl","shell.ampl"]
        if(psf_gauss.find("dc") != -1):
            param_names          = ["r0","xpos","ypos","xrays_mod.ampl","disc.ampl"]
    else:
        if (psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            param_names = ["fwhm", "xpos", "ypos", "ellip", "theta"]
        if (psf_gauss.find("2G") != -1):
            param_names = ["fwhm", "xpos", "ypos"]
        if (psf_gauss.find("G3") != -1 or psf_gauss.find("G2") != -1):
            param_names = ["fwhm", "xpos", "ypos"]
        if (psf_gauss.find("sh") != -1):
            param_names = ["r0", "xpos", "ypos", "width"]
        if (psf_gauss.find("dc") != -1):
            param_names = ["r0", "xpos", "ypos"]

    data                 = []
    params               = []
    val                  = np.zeros((dof,3))
    val_min              = np.zeros((dof,3))
    val_max              = np.zeros((dof,3))

    f                    = open(saved_results)
    data                 = np.ndarray((len(column),dof))
    
    for l in range (0,len(conf_names)):
        temp   = f.readline().split()
        if (temp[0].find(column[0]) != -1):
                params.append(temp[2].split("(")[1].split(",")[0])
                for j in range (3,len(temp)-1):
                    params.append(temp[j].split()[0].split(",")[0])
                params.append(temp[len(temp)-1].split(")")[0])
        for k in range(1,len(column)):
            if (temp[0].find(column[k]) != -1):
                for j in range(2,len(temp)):
                    if (j == 2):
                        if( temp[j].split("(")[1].split(",")[0] == None):
                            temp[j].split("(")[1].split(",")[0] = 0
                        data[k][j-2] = float(temp[j].split("(")[1].split(",")[0])
                    if (j == len(temp)-1):
                        if( temp[j].split(")")[0] == None):
                            temp[j].split(")")[0] = 0
                        data[k][j-2] = float(temp[j].split(")")[0])
                    if (j > 2 and j < len(temp)-1):
                        if (temp[j].split()[0].split(",")[0] == None):
                            temp[j].split()[0].split(",")[0] = 0
                        data[k][j-2] = float(temp[j].split()[0].split(",")[0])
                                            
    print("\n Parameter results (from the image):\n")
    for i in range(0,dof):
        print "Parameter :", params[i]," = ",data[1][i],"  ",data[2][i]," + ",data[3][i],"\n"


    for j in range(0,len(param_names)):
        k = 0
        for i in range(1,len(params)):
            if(params[i].find(param_names[j]) != -1):
                val[j][k]     = (data[1][i])
                val_min[j][k] = (data[2][i])
                val_max[j][k] = (data[3][i])
                k             = k + 1


    if (save_fit.find("X") != -1):

        if (psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            sys.stdout.write("\n Coordinates of the asymetrical gaussian :\n")
            alpha_gauss1, delta_gauss1 = conversion.convert_pix2wcs(fits_filename, val[1][0], val[2][0])
            conversion.convert_deg2timearcs(alpha_gauss1, delta_gauss1)
            print "Error for the RA  : ", conversion.error_convert_ra(fits_filename,
                                                                      val_min[1][0]), conversion.error_convert_ra(
                fits_filename, val_max[1][0])
            print "Error for the DEC : ", conversion.error_convert_dec(fits_filename,
                                                                       val_min[2][0]), conversion.error_convert_dec(
                fits_filename, val_min[2][0])

            sys.stdout.write("\n Source morphology \n")
            print "\n Intrinsic standart deviation (major axis)  : ", conversion.fwhm_sigma(val[0][0],
                                                                                            factor_pix2deg_gamma), "  ", conversion.fwhm_sigma(
                val_min[0][0], factor_pix2deg_gamma), " + ", conversion.fwhm_sigma(val_max[0][0],
                                                                                   factor_pix2deg_gamma), " arcmin"
            print "\n Intrinsic standart deviation (minor axis)  : ", conversion.fwhm_sigma(val[0][0] * (1 - val[3][0]),
                                                                                            factor_pix2deg_gamma), " arcmin"
            print "\n Inclination                                :", val[4][0] * 180. / np.pi, "  ", val_min[4][
                                                                                                         0] * 180. / np.pi, " + ", \
            val_max[4][0] * 180. / np.pi, " degrees"

            print "\nComponent amplitude / X-ray amplitude :", val[6][0] / val[5][0]

        if (psf_gauss.find("sh") != -1):
            sys.stdout.write("\n Coordinates of the shell :\n")
            alpha_gauss1, delta_gauss1 = conversion.convert_pix2wcs(fits_filename, val[1][0], val[2][0])
            conversion.convert_deg2timearcs(alpha_gauss1, delta_gauss1)
            print "Error for the RA  : ", conversion.error_convert_ra(fits_filename,
                                                                      val_min[1][0]), conversion.error_convert_ra(
                fits_filename, val_max[1][0])
            print "Error for the DEC : ", conversion.error_convert_dec(fits_filename,
                                                                       val_min[2][0]), conversion.error_convert_dec(
                fits_filename, val_min[2][0])

            sys.stdout.write("\n Shell morphology \n")
            print "\n R      : ", val[0][0] * factor_pix2deg_gamma * 60., " - ", val_min[0][
                                                                                     0] * factor_pix2deg_gamma * 60., " + ", \
            val_max[0][0] * factor_pix2deg_gamma * 60., " arcmin"
            print "\n Width  : ", val[3][0] * factor_pix2deg_gamma * 60., " - ", val_min[3][
                                                                                     0] * factor_pix2deg_gamma * 60., " + ", \
            val_max[3][0] * factor_pix2deg_gamma * 60., " arcmin"

            print "\nComponent amplitude / X-ray amplitude :", val[5][0] / val[4][0]

        if (psf_gauss.find("dc") != -1):
            sys.stdout.write("\n Coordinates of the disk :\n")
            alpha_gauss1, delta_gauss1 = conversion.convert_pix2wcs(fits_filename, val[1][0], val[2][0])
            conversion.convert_deg2timearcs(alpha_gauss1, delta_gauss1)
            print "Error for the RA  : ", conversion.error_convert_ra(fits_filename,
                                                                      val_min[1][0]), conversion.error_convert_ra(
                fits_filename, val_max[1][0])
            print "Error for the DEC : ", conversion.error_convert_dec(fits_filename,
                                                                       val_min[2][0]), conversion.error_convert_dec(
                fits_filename, val_min[2][0])

            sys.stdout.write("\n Disk parameters \n")
            print "\n R      : ", val[0][0] * factor_pix2deg_gamma * 60., " - ", val_min[0][
                                                                                     0] * factor_pix2deg_gamma * 60., " + ", \
            val_max[0][0] * factor_pix2deg_gamma * 60., " arcmin"

            print "\nComponent amplitude / X-ray amplitude :", val[4][0] / val[3][0]


        if (psf_gauss.find("2G") != -1):
            sys.stdout.write("\n Coordinates of the symetrical gaussian :\n")
            alpha_gauss1, delta_gauss1 = conversion.convert_pix2wcs(fits_filename, val[1][0], val[2][0])
            conversion.convert_deg2timearcs(alpha_gauss1, delta_gauss1)
            print "Error for the RA  : ", conversion.error_convert_ra(fits_filename,
                                                                      val_min[1][0]), conversion.error_convert_ra(
                fits_filename, val_max[1][0])
            print "Error for the DEC : ", conversion.error_convert_dec(fits_filename,
                                                                       val_min[2][0]), conversion.error_convert_dec(
                fits_filename, val_min[2][0])

            sys.stdout.write("\n Source morphology \n")
            print "\n Intrinsic standart deviation   : ", conversion.fwhm_sigma(val[0][0],
                                                                                factor_pix2deg_gamma), "  ", conversion.fwhm_sigma(
                val_min[0][0], factor_pix2deg_gamma), " + ", conversion.fwhm_sigma(val_max[0][0],
                                                                                   factor_pix2deg_gamma), " arcmin"
            print "\nComponent amplitude / X-ray amplitude :", val[4][0] / val[3][0]

    else:

        if(psf_gauss.find("G3") != -1):
            sys.stdout.write("\n Coordinates of the gaussian on the left :\n")
            alpha_gauss1,delta_gauss1   = conversion.convert_pix2wcs(fits_filename, val[1][0], val[2][0])
            conversion.convert_deg2timearcs(alpha_gauss1,delta_gauss1)
            print "Error for the RA  : ",conversion.error_convert_ra(fits_filename,val_min[1][0]),conversion.error_convert_ra(fits_filename,val_max[1][0])
            print "Error for the DEC : ",conversion.error_convert_dec(fits_filename,val_min[2][0]),conversion.error_convert_dec(fits_filename,val_min[2][0])

        
            sys.stdout.write("\n Coordinates of the gaussian on the right :\n")
            alpha_gauss2,delta_gauss2   = conversion.convert_pix2wcs(fits_filename, val[1][1], val[2][1])
            conversion.convert_deg2timearcs(alpha_gauss2,delta_gauss2)
            print "Error for the RA  : ",conversion.error_convert_ra(fits_filename,val_min[1][1]),conversion.error_convert_ra(fits_filename,val_max[1][1])
            print "Error for the DEC : ",conversion.error_convert_dec(fits_filename,val_min[2][1]),conversion.error_convert_dec(fits_filename,val_min[2][1])

            sys.stdout.write("\n Coordinates of the gaussian at the center :\n")
            alpha_gauss3,delta_gauss3   = conversion.convert_pix2wcs(fits_filename, val[1][2], val[2][2])
            conversion.convert_deg2timearcs(alpha_gauss3,delta_gauss3)
            print "Error for the RA  : ",conversion.error_convert_ra(fits_filename,val_min[1][2]),conversion.error_convert_ra(fits_filename,val_max[1][2])
            print "Error for the DEC : ",conversion.error_convert_dec(fits_filename,val_min[2][2]),conversion.error_convert_dec(fits_filename,val_min[2][2])
        
            sys.stdout.write("\n Source morphology \n")
        
            print "\n Intrinsic standart deviation (left gaussian)  : ",conversion.fwhm_sigma(val[0][0],factor_pix2deg_gamma)," - ",conversion.fwhm_sigma(val_min[0][0],factor_pix2deg_gamma)," + ",conversion.fwhm_sigma(val_max[0][0],factor_pix2deg_gamma)," arcmin"
            print "\n Intrinsic standart deviation (right gaussian) : ",conversion.fwhm_sigma(val[0][1],factor_pix2deg_gamma)," - ",conversion.fwhm_sigma(val_min[0][1],factor_pix2deg_gamma)," + ",conversion.fwhm_sigma(val_max[0][1],factor_pix2deg_gamma)," arcmin"
            print "\n Intrinsic standart deviation (center)         : ",conversion.fwhm_sigma(val[0][2],factor_pix2deg_gamma)," - ",conversion.fwhm_sigma(val_min[0][2],factor_pix2deg_gamma)," + ",conversion.fwhm_sigma(val_max[0][2],factor_pix2deg_gamma)," arcmin"

        if(psf_gauss.find("G2") != -1):
            sys.stdout.write("\n Coordinates of the gaussian on the left :\n")
            alpha_gauss1,delta_gauss1   = conversion.convert_pix2wcs(fits_filename, val[1][0], val[2][0])
            conversion.convert_deg2timearcs(alpha_gauss1,delta_gauss1)
            print "Error for the RA  : ",conversion.error_convert_ra(fits_filename,val_min[1][0]),conversion.error_convert_ra(fits_filename,val_max[1][0])
            print "Error for the DEC : ",conversion.error_convert_dec(fits_filename,val_min[2][0]),conversion.error_convert_dec(fits_filename,val_min[2][0])
    
            sys.stdout.write("\n Coordinates of the gaussian on the right :\n")
            alpha_gauss2,delta_gauss2   = conversion.convert_pix2wcs(fits_filename, val[1][1], val[2][1])
            conversion.convert_deg2timearcs(alpha_gauss2,delta_gauss2)
            print "Error for the RA  : ",conversion.error_convert_ra(fits_filename,val_min[1][1]),conversion.error_convert_ra(fits_filename,val_max[1][1])
            print "Error for the DEC : ",conversion.error_convert_dec(fits_filename,val_min[2][1]),conversion.error_convert_dec(fits_filename,val_min[2][1])
        
            sys.stdout.write("\n Source morphology \n")
            print "\n Intrinsic standart deviation (left gaussian)  : ",conversion.fwhm_sigma(val[0][0],factor_pix2deg_gamma)," - ",conversion.fwhm_sigma(val_min[0][0],factor_pix2deg_gamma)," + ",conversion.fwhm_sigma(val_max[0][0],factor_pix2deg_gamma)," arcmin"
            print "\n Intrinsic standart deviation (right gaussian) : ",conversion.fwhm_sigma(val[0][1],factor_pix2deg_gamma)," - ",conversion.fwhm_sigma(val_min[0][1],factor_pix2deg_gamma)," + ",conversion.fwhm_sigma(val_max[0][1],factor_pix2deg_gamma)," arcmin"
    

    if(nb_displays == 2):
        sys.stdout.close()
        sys.stdout = temporary   
        
    return 


def quick_results(save_fit, psf_gauss, name_file, elliptical,alpha,stat,dof):

    if (save_fit.find("yes") != -1):
        nb_displays = 2
        temporary   = sys.stdout

        sys.stdout  = open(names.results_path + "quick_fit_results" + psf_gauss + ".txt", 'w')
    else :
        nb_displays = 1
        temporary   = sys.stdout


    if(names.analysis == 0):
        analysis = 'hgps'
    else:
        analysis = 'pa'

    ra_psr   = 228.5292
    dec_psr  = -59.1575
    glon_psr = 320.34
    glat_psr = -1.20


    print "alpha   : ",alpha,"\n"
    print "CASH    : ",stat,"\n"
    print "d.o.f   : ",dof + 1,"\n"
    print "AIC     : ",stat + 2 * (dof + 1),"\n"


    if(psf_gauss.find("G") != -1):
        xray_ampl  = elliptical[5]
        comp_ampl  = elliptical[6]
        comp_xpos  = elliptical[0]
        comp_ypos  = elliptical[1]
        comp_ellip = elliptical[3]
        comp_theta = elliptical[4]
        comp_fwhm  = elliptical[2]


        print "X           : ",comp_xpos
        print "Y           : ",comp_ypos
        print "ellipticity : ",comp_ellip
        print "theta       : ",comp_theta
        print "FWHM        : ",comp_fwhm,"\n"


    if(psf_gauss.find("sh") != -1):
        xray_ampl  = elliptical[5]
        comp_ampl  = elliptical[4]
        comp_xpos  = elliptical[0]
        comp_ypos  = elliptical[1]
        comp_r0    = elliptical[2]
        comp_width = elliptical[3]

        print "X           : ",comp_xpos
        print "Y           : ",comp_ypos
        print "r0          : ",comp_r0
        print "width       : ",comp_width,"\n"


    if(psf_gauss.find("dc") != -1):
        xray_ampl  = elliptical[5]
        comp_ampl  = elliptical[3]
        comp_xpos  = elliptical[0]
        comp_ypos  = elliptical[1]
        comp_r0    = elliptical[2]

        print "X           : ",comp_xpos
        print "Y           : ",comp_ypos
        print "r0          : ",comp_r0,"\n"


    if(psf_gauss.find("Sc") != -1):
        xray_ampl  = elliptical[5]
        comp_ampl  = elliptical[6]
        comp_xpos  = elliptical[0]
        comp_ypos  = elliptical[1]
        comp_ellip = elliptical[3]
        comp_theta = elliptical[4]
        comp_fwhm  = elliptical[2]
        comp_n     = elliptical[7]

        print "X           : ",comp_xpos
        print "Y           : ",comp_ypos
        print "ellipticity : ",comp_ellip
        print "theta       : ",comp_theta
        print "r0          : ",comp_fwhm
        print "Sersic n    : ",comp_n,"\n"

    if(psf_gauss.find("G") == -1 and psf_gauss.find("sh") == -1 and psf_gauss.find("dc") == -1 and psf_gauss.find("Sc") == -1):
        xray_ampl = elliptical[0]
        print "X-ray amplitude : ", xray_ampl

    else:
        print "Component amplitude / X-ray amplitude : ", comp_ampl / xray_ampl

        if (analysis.find("hgps") != -1):
            binsize    = 0.02
            pic_center = 100.5

            coord1, coord2 = conversion.convert_image(glon_psr, glat_psr, pic_center, comp_xpos, comp_ypos, binsize)

            equ = SkyCoord(l=coord1, b=coord2, frame='galactic', unit='deg').fk5

            print "l : ", coord1, " , b : ", coord2


        if (analysis.find("pa") != -1):
            binsize    = 0.01
            pic_center = 200.5

            coord1, coord2 = conversion.convert_image(ra_psr, dec_psr, pic_center, comp_xpos, comp_ypos, binsize)

            equ = SkyCoord(ra=coord1, dec=coord2, unit='deg')

            print "l : ", equ.galactic.l.deg, " , b : ", equ.galactic.b.deg

        conversion.convert_deg2timearcs(equ.ra.deg, equ.dec.deg)

        if (psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):

            print "Major axis :", conversion.fwhm_sigma(comp_fwhm, binsize),  "\n"

            print "Minor axis :", conversion.fwhm_sigma(comp_fwhm * (1 - comp_ellip), binsize), "\n"

            if (comp_theta >= 0.):
                print "Inclination : ", comp_theta * 180. / np.pi, "\n"
            if (comp_theta < 0.):
                print "Inclination : ", (360. + (comp_theta * 180. / np.pi)), "\n"


        if (psf_gauss.find("2G") != -1):

            print "Intrinsic s :", conversion.fwhm_sigma(comp_fwhm, binsize),  "\n"


        if(psf_gauss.find("sh") != -1):

            print "Inner radius : ", conversion.fwhm_sigma(comp_r0, binsize),  "\n"
            print "Outer radius : ", conversion.fwhm_sigma(comp_r0 + comp_width, binsize),  "\n"
            print "Width        : ", conversion.fwhm_sigma(comp_width, binsize),  "\n"


        if(psf_gauss.find("dc") != -1):
            print "Disk radius  : ", conversion.fwhm_sigma(comp_r0, binsize),  "\n"


    if(nb_displays == 2):
        sys.stdout.close()
        sys.stdout = temporary


    return


def quick_errors(save_fit_true,save_fit, psf_gauss,elliptical, configuration):
    if (save_fit_true.find("yes") != -1):
        nb_displays = 2
        temporary = sys.stdout
        sys.stdout = open(names.results_path + "quick_errors_" + psf_gauss + ".txt", 'w')
    else:
        nb_displays = 1
        temporary = sys.stdout

    ra_psr   = 228.5292
    dec_psr  = -59.1575
    glon_psr = 320.34
    glat_psr = -1.20

    if (save_fit.find("xr") != -1):
        saved_results = names.save_path[names.analysis] + "conf_" + psf_gauss + ".txt"
        saved_fitinfo = names.save_path[names.analysis] + "fit_info_" + psf_gauss[:3] + ".out"


    print("\n ***********\n **********\n ________________________\n ####### Results ######\n")
    print("Object      : " + names.source[names.analysis] + "\n")
    print("Fit : " + names.tags[names.analysis] + " with " + names.fit_component[configuration])
    print("Alpha parameter : " + str(names.alpha_factor_pa[configuration]))

    nf, cash_stats = np.loadtxt(saved_fitinfo, dtype=float, usecols=(0, 1), unpack=True, skiprows=1)
    cash = cash_stats[-1]

    with file(saved_fitinfo) as f:
        one_row = f.readline()
    dof = len(one_row.split()) - 3
    f.close()

    with open(saved_results, 'r') as r_file:
        filedata = r_file.read()
    filedata = filedata.replace('None', '0')

    with open(saved_results, 'w') as r_file:
        r_file.write(filedata)
    r_file.close()

    if (save_fit.find("xr0") == -1):
        print "\n STATISTICS :\n   CASH  = ", cash, "\n   d.o.f : ", dof ,"+ 1", "\n   AIC   = ", cash + 2. * (dof + 1)
    if(save_fit.find("xr0") != -1):
        print "\n STATISTICS :\n   CASH  = ", cash, "\n   d.o.f : ", dof, "\n   AIC   = ", cash + 2. * dof
    if (save_fit.find("xr") == -1):
        print "\n STATISTICS :\n   CASH  = ", cash, "\n   d.o.f : ", dof, "\n   AIC   = ", cash + 2. * dof


    conf_names, equals = np.loadtxt(saved_results, dtype=str, usecols=(0, 1), unpack=True, skiprows=1)
    column = ["parnames", "parvals", "parmins", "parmaxes"]

    if (save_fit.find("x") != -1):
        if (psf_gauss.find('G') == -1 and psf_gauss.find('sh') == -1 and psf_gauss.find('dc') == -1):
            param_names = ["ampl"]

        if (psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            param_names = ["fwhm", "xpos", "ypos", "ellip", "theta", "xrays_mod.ampl", "gcomp.ampl"]
        if (psf_gauss.find("2G") != -1):
            param_names = ["fwhm", "xpos", "ypos", "xrays_mod.ampl", "gcomp.ampl"]
        if (psf_gauss.find("sh") != -1):
            param_names = ["r0", "xpos", "ypos", "width", "xrays_mod.ampl", "shell.ampl"]
        if (psf_gauss.find("dc") != -1):
            param_names = ["r0", "xpos", "ypos", "xrays_mod.ampl", "disc.ampl"]
        if (psf_gauss.find("Sc") != -1 and psf_gauss.find("P") == -1):
            param_names = ["r0", "xpos", "ypos", "ellip", "theta", "xrays_mod.ampl", "gcomp.ampl", "gcomp.n"]
        if (psf_gauss.find("Sc") != -1 and psf_gauss.find("P") != -1):
            param_names = ["r0", "ellip", "theta", "xrays_mod.ampl", "gcomp.ampl", "gcomp.n"]

    else:
        if (psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            param_names = ["fwhm", "xpos", "ypos", "ellip", "theta"]
        if (psf_gauss.find("2G") != -1):
            param_names = ["fwhm", "xpos", "ypos"]
        if (psf_gauss.find("G3") != -1 or psf_gauss.find("G2") != -1):
            param_names = ["fwhm", "xpos", "ypos"]
        if (psf_gauss.find("sh") != -1):
            param_names = ["r0", "xpos", "ypos", "width"]
        if (psf_gauss.find("dc") != -1):
            param_names = ["r0", "xpos", "ypos"]

    data    = []
    params  = []
    val     = np.zeros((dof, 3))
    val_min = np.zeros((dof, 3))
    val_max = np.zeros((dof, 3))

    f = open(saved_results)
    data = np.ndarray((len(column), dof))

    for l in range(0, len(conf_names)):
        temp = f.readline().split()
        if (temp[0].find(column[0]) != -1):
            params.append(temp[2].split("(")[1].split(",")[0])
            for j in range(3, len(temp) - 1):
                params.append(temp[j].split()[0].split(",")[0])
            params.append(temp[len(temp) - 1].split(")")[0])
        for k in range(1, len(column)):
            if (temp[0].find(column[k]) != -1):
                for j in range(2, len(temp)):
                    if (j == 2):
                        if (temp[j].split("(")[1].split(",")[0].find('None') != -1 or temp[j].split("(")[1].split(",")[0] == None or type(temp[j].split(",")[0]) == str):
                            temp[j].split("(")[1].split(",")[0] = 0
                        data[k][j - 2] = float(temp[j].split("(")[1].split(",")[0])
                    if (j == len(temp) - 1 and j != 2):
                        if (temp[j].split(")")[0].find('None') != -1 or temp[j].split(")")[0] == None or type(temp[j].split(")")) == str):
                            temp[j].split(")")[0] = 0
                        data[k][j - 2] = float(temp[j].split(")")[0])
                    if (j > 2 and j < len(temp) - 1):
                        if (temp[j].split()[0].split(",")[0].find('None') != -1 or temp[j].split()[0].split(",")[0] == None or type(temp[j].split()[0].split(",")[0]) == str):
                            temp[j].split()[0].split(",")[0] = 0

                        data[k][j - 2] = float(temp[j].split()[0].split(",")[0])

    print("\n Parameter results (from the image):\n")
    for i in range(0, dof):
        print "Parameter :", params[i], " = ", data[1][i], "  ", data[2][i], " + ", data[3][i], "\n"

    if( len(param_names) != 1):
        for j in range(0, len(param_names)):
            k = 0
            for i in range(0, len(params)):
                if (params[i].find(param_names[j]) != -1):
                    val[j][k] = (data[1][i])
                    val_min[j][k] = (data[2][i])
                    val_max[j][k] = (data[3][i])
                    k = k + 1
    else :
        for j in range(0, len(param_names)):
            k = 0
            for i in range(0, 1):
                if(params[i].find(param_names[j]) != -1):

                    val[j][k]     = (data[1][i])
                    val_min[j][k] = (data[2][i])
                    val_max[j][k] = (data[3][i])



    if (save_fit.find("x") != -1):

        if (psf_gauss.find("G") != -1):
            xray_ampl  = elliptical[5]
            comp_ampl  = elliptical[6]
            comp_xpos  = elliptical[0]
            comp_ypos  = elliptical[1]
            comp_ellip = elliptical[3]
            comp_theta = elliptical[4]
            comp_fwhm  = elliptical[2]

        if (psf_gauss.find("sh") != -1):
            xray_ampl  = elliptical[5]
            comp_ampl  = elliptical[4]
            comp_xpos  = elliptical[0]
            comp_ypos  = elliptical[1]
            comp_r0    = elliptical[2]
            comp_width = elliptical[3]


        if (psf_gauss.find("dc") != -1):
            xray_ampl = elliptical[5]
            comp_ampl = elliptical[3]
            comp_xpos = elliptical[0]
            comp_ypos = elliptical[1]
            comp_r0   = elliptical[2]

        if (psf_gauss.find("Sc") != -1):
            xray_ampl = elliptical[5]
            comp_ampl = elliptical[6]
            comp_xpos = elliptical[0]
            comp_ypos = elliptical[1]
            comp_ellip = elliptical[3]
            comp_theta = elliptical[4]
            comp_fwhm  = elliptical[2]
            comp_n     = elliptical[7]

        if (psf_gauss.find("G") == -1 and psf_gauss.find("sh") == -1 and psf_gauss.find("dc") == -1 and psf_gauss.find("Sc") == -1):
            xray_ampl = elliptical[0]

        else:

            if (psf_gauss.find("Sc") != -1 and psf_gauss.find("P") != -1):
                binsize    = 0.01
                pic_center = 200.5

                coord1, coord2     = conversion.convert_image(ra_psr, dec_psr, pic_center, comp_xpos, comp_ypos, binsize)

                equ            = SkyCoord(ra=coord1, dec=coord2, unit='deg')

                print "Frozen PSR coordinates."

                print "l : ", equ.galactic.l.deg
                print "b : ", equ.galactic.b.deg

                conversion.convert_deg2timearcs(equ.ra.deg, equ.dec.deg)

                print "X-ray amplitude ", xray_ampl, "- / + ", val_min[0][0], val_max[0][0]

                print "\nComponent amplitude / X-ray amplitude :", val[4][0] / val[3][0]," - / + ",val[4][0] / val[3][0] * (val_min[4][0] / val[3][0] + val_min[3][0]/val[3][0]),val[4][0] / val[3][0] * (val_max[4][0] / val[4][0] + val_max[3][0]/val[3][0]),"\n"

                sys.stdout.write("\n Characteristics of the Sersic profile :\n")

                print "Sersic index n = ",comp_n

                print "Major axis :", conversion.fwhm_sigma(comp_fwhm, binsize), " - / +",conversion.fwhm_sigma(val_min[0][0],binsize),conversion.fwhm_sigma(val_max[0][0],binsize),"arcmin\n"

                print "Minor axis :", conversion.fwhm_sigma(comp_fwhm * (1 - comp_ellip), binsize),"arcmin\n"

                if (comp_theta >= 0.):
                    print "Inclination : ", comp_theta * 180. / np.pi,"- / + ", val_min[2][0] * 180. / np.pi,val_max[2][0] * 180. / np.pi,"degrees \n"
                if (comp_theta < 0.):
                    print "Inclination : ", (360. + (comp_theta * 180. / np.pi)),"- / + ",(360. + (val_min[2][0] * 180. / np.pi)),(360. + (val_max[2][0] * 180. / np.pi)), "degrees\n"

                print "Sersic profile : "
                sersic   = conversion.Sersic_profile(comp_n, comp_fwhm, comp_ellip, comp_theta, comp_xpos, comp_ypos, comp_ampl)


            else :

                binsize    = 0.01
                pic_center = 200.5

                coord1, coord2     = conversion.convert_image(ra_psr, dec_psr, pic_center, comp_xpos, comp_ypos, binsize)

                ref1, ref2         = conversion.convert_image(ra_psr,dec_psr,pic_center,0.,0.,binsize)
                err1_min, err2_min = conversion.convert_image(ra_psr, dec_psr, pic_center, val_min[1][0], 0., binsize)
                error1_min         = abs(err1_min - ref1)
                err1_min, err2_min = conversion.convert_image(ra_psr, dec_psr, pic_center, 0., val_min[2][0], binsize)
                error2_min         = abs(err2_min - ref2)

                err1_max, err2_max = conversion.convert_image(ra_psr, dec_psr, pic_center, val_max[1][0], 0.,binsize)
                error1_max         = abs(err1_max - ref1)
                err1_max, err2_max = conversion.convert_image(ra_psr, dec_psr, pic_center, 0., val_max[2][0],binsize)
                error2_max         = abs(err2_max - ref2)

                equ            = SkyCoord(ra=coord1, dec=coord2, unit='deg')
                equ_err_min    = SkyCoord(ra=error1_min, dec=error2_min, unit='deg')
                equ_err_max    = SkyCoord(ra=error1_max, dec=error2_max, unit='deg')

                print "l : ", equ.galactic.l.deg," - ",equ_err_min.galactic.l.deg," + ",equ_err_max.galactic.l.deg
                print "b : ", equ.galactic.b.deg," - ",equ_err_min.galactic.b.deg," + ",equ_err_max.galactic.b.deg

                conversion.convert_deg2timearcs(equ.ra.deg, equ.dec.deg)
                print "Error for the RA  : ", equ_err_min.ra.deg *12*3600/180., equ_err_max.ra.deg * 12 * 3600 / 180.
                print "Error for the DEC : ", equ_err_min.dec.deg *3600, equ_err_max.dec.deg * 3600

                print "X-ray amplitude ", xray_ampl, "- / + ", val_min[0][0], val_max[0][0]

                if (psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):

                    print "\nComponent amplitude / X-ray amplitude :", val[6][0] / val[5][0]," - / + ",val[6][0] / val[5][0] * (val_min[6][0] / val[6][0] + val_min[5][0]/val[5][0]),val[6][0] / val[5][0] * (val_max[6][0] / val[6][0] + val_max[5][0]/val[5][0]),"\n"

                    sys.stdout.write("\n Coordinates of the asymetrical gaussian :\n")

                    print "Major axis :", conversion.fwhm_sigma(comp_fwhm, binsize), " - / +",conversion.fwhm_sigma(val_min[0][0],binsize),conversion.fwhm_sigma(val_max[0][0],binsize),"arcmin\n"

                    print "Minor axis :", conversion.fwhm_sigma(comp_fwhm * (1 - comp_ellip), binsize),"arcmin\n"

                    if (comp_theta >= 0.):
                        print "Inclination : ", comp_theta * 180. / np.pi,"- / + ", val_min[4][0] * 180. / np.pi,val_max[4][0] * 180. / np.pi,"degrees \n"
                    if (comp_theta < 0.):
                        print "Inclination : ", (360. + (comp_theta * 180. / np.pi)),"- / + ",(360. + (val_min[4][0] * 180. / np.pi)),(360. + (val_max[4][0] * 180. / np.pi)), "degrees\n"

                if (psf_gauss.find("2G") != -1):
                    # Carefull!! the error bars for the fraction are not well computated
                    print "\nComponent amplitude / X-ray amplitude :", val[4][0] / (val[3][0] + val[4][0]),\
                        " - / + ",val[4][0] / val[3][0] * (val_min[4][0] / val[4][0] + val_min[3][0]/val[3][0]),val[4][0] / val[3][0] * (val_max[4][0] / val[4][0] + val_max[3][0]/val[3][0]),"\n"

                    intG   = val[4][0] * 2. * np.pi * (comp_fwhm / (2. * np.sqrt(2. * np.log(2.))))**2
                    intG_X = intG/ ( val[3][0] + intG)
                    ern_A_ = val_min[4][0] / val[4][0]
                    ern_X  = val_min[3][0] / val[3][0]
                    ern_s2 = 2. * (comp_fwhm / (2. * np.sqrt(2. * np.log(2.))) / (val_min[0][0])/ (2. * np.sqrt(2. * np.log(2.))))
                    erp_A  = val_max[4][0] / val[4][0]
                    erp_X  = val_max[3][0] / val[3][0]
                    erp_s2 = 2. * (comp_fwhm / (2. * np.sqrt(2. * np.log(2.))) / (val_max[0][0])/ (2. * np.sqrt(2. * np.log(2.))))

                    print "\nIntegral of the Gaussian / X-ray amplitude :   ", intG_X,\
                        " - / + ",intG_X * ( ern_A_+ ern_X + ern_s2)**-1 ,\
                        intG_X * (erp_A + erp_X + erp_s2)**-1,"\n"

                    print "Intrinsic s :", conversion.fwhm_sigma(comp_fwhm, binsize), " - / +",conversion.fwhm_sigma(val_min[0][0],binsize),conversion.fwhm_sigma(val_max[0][0],binsize), "\n"

                if (psf_gauss.find("sh") != -1):

                    print "\nComponent amplitude / X-ray amplitude :", val[5][0] / val[4][0]," - / + ",val[5][0] / val[4][0] * (val_min[5][0] / val[5][0] + val_min[4][0]/val[4][0]),val[5][0] / val[4][0] * (val_max[5][0] / val[5][0] + val_max[4][0]/val[4][0]),"\n"

                    print "Inner radius : ", conversion.fwhm_sigma(comp_r0, binsize),  " - / +",conversion.fwhm_sigma(val_min[0][0],binsize),conversion.fwhm_sigma(val_max[0][0],binsize),"\n"
                    print "Outer radius : ", conversion.fwhm_sigma(comp_r0 + comp_width, binsize),  " - / +",conversion.fwhm_sigma(val_min[0][0] + val_min[3][0],binsize),conversion.fwhm_sigma(val_max[0][0] + val_max[3][0],binsize),"\n"
                    print "Width        : ", conversion.fwhm_sigma(comp_width, binsize),  " - / +",conversion.fwhm_sigma(val_min[3][0],binsize),conversion.fwhm_sigma(val_max[3][0],binsize),"\n"

                if (psf_gauss.find("dc") != -1):

                    print "\nComponent amplitude / X-ray amplitude :", val[4][0] / val[3][0]," - / + ",val[4][0] / val[3][0] * (val_min[4][0] / val[4][0] + val_min[3][0]/val[3][0]),val[4][0] / val[3][0] * (val_max[4][0] / val[4][0] + val_max[3][0]/val[3][0]),"\n"

                    print "Disk radius  : ", conversion.fwhm_sigma(comp_r0, binsize), " - / +",conversion.fwhm_sigma(val_min[0][0],binsize),conversion.fwhm_sigma(val_max[0][0],binsize),"arcmin\n"


                if (psf_gauss.find("Sc") != -1):

                    print "\nComponent amplitude / X-ray amplitude :", val[6][0] / val[5][0]," - / + ",val[6][0] / val[5][0] * (val_min[6][0] / val[6][0] + val_min[5][0]/val[5][0]),val[6][0] / val[5][0] * (val_max[6][0] / val[6][0] + val_max[5][0]/val[5][0]),"\n"

                    sys.stdout.write("\n Coordinates of the Sersic component :\n")

                    print "Sersic index n : ", val[7][0]," - / +",val_min[7][0],val_max[7][0]

                    print "Major axis :", conversion.fwhm_sigma(comp_fwhm, binsize), " - / +",conversion.fwhm_sigma(val_min[0][0],binsize),conversion.fwhm_sigma(val_max[0][0],binsize),"arcmin\n"

                    print "Minor axis :", conversion.fwhm_sigma(comp_fwhm * (1 - comp_ellip), binsize),"arcmin\n"

                    if (comp_theta >= 0.):
                        print "Inclination : ", comp_theta * 180. / np.pi,"- / + ", val_min[4][0] * 180. / np.pi,val_max[4][0] * 180. / np.pi,"degrees \n"
                    if (comp_theta < 0.):
                        print "Inclination : ", (360. + (comp_theta * 180. / np.pi)),"- / + ",(360. + (val_min[4][0] * 180. / np.pi)),(360. + (val_max[4][0] * 180. / np.pi)), "degrees\n"

                    print "Sersic profile : "
                    sersic   = conversion.Sersic_profile(comp_n, comp_fwhm, comp_ellip, comp_theta, comp_xpos, comp_ypos, comp_ampl)

        if (nb_displays == 2):
            sys.stdout.close()
            sys.stdout = temporary

    return

