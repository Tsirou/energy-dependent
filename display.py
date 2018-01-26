import os
import names

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
import scipy.ndimage as ndimage


from astropy import wcs
from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D

from names import path,filename_gamma
from matplotlib.patches import Ellipse,Circle

from conversion import fwhm_sigma,affine

#==============================================================================
# For cool display stuff
#==============================================================================




#==============================================================================
#3D scatter plot
def plot_3d_lin(xx,yy,zz,ax,variable_names,colors,names,title):
    ax.set_figsize=(30000,30000) 
    ax.set_title(title,fontsize=30,y=1.08)  

    ax.set_xlabel(variable_names[0],fontsize=25,y=1.7)    
    ax.set_ylabel(variable_names[1],fontsize=25,x=1.7)
    ax.set_zlabel(variable_names[2],fontsize=25,y=1.7)
    
    ax.tick_params(axis='both',which='major',labelsize=20) 

    scat = ax.scatter(xx,yy,zz, zdir='z',c=colors,alpha=0.5,label=names)
        
    return scat
        
# 2D with colorbar with a 3D projection
def plot_2d_cmap(xx,yy,zz,fig3,variable_names,title):
    ax3 = fig3.add_subplot(111,projection='3d')

    ax3.set_figsize=(3000,3000)    
    
    ax3.set_title(title,fontsize=60,y=1.08)    
    
    ax3.set_xlabel(variable_names[0],fontsize=25)
    ax3.set_ylabel(variable_names[1],fontsize=25)

    p  = ax3.scatter(xx,yy,c = zz,zdir='z',cmap = "hot")

    ax3.w_zaxis.line.set_lw(0.)
    ax3.set_zticks([])
    ax3.view_init(azim=180,elev=-90)
    fig3.colorbar(p,label=variable_names[2])
    

# "Chosen" features display
def plot_display(x,y,dpis):
    #Size of the figure
    fig = plt.figure(figsize=(x,y), dpi=dpis)
    mpl.rc('font',family='Serif')
    mpl.rcParams.update({'font.size': 25})
    mpl.rcParams.update({'xtick.labelsize': 30})
    mpl.rcParams.update({'ytick.labelsize': 30})
    
    return fig

# Reset to default display configuration
def reset_display():
    mpl.rcParams.update(mpl.rcParamsDefault)

# Refreshing the scatter plot
def update_plot(i, data, scat):

    scat.set_array(data[i])

    return scat

# JPEG file generation [From the web]
def make_views(ax,angles,elevation=None, width=20, height = 20,
                prefix='tmprot_',**kwargs):
    """
    Makes jpeg pictures of the given 3d ax, with different angles.
    Args:
        ax (3D axis): te ax
        angles (list): the list of angles (in degree) under which to
                       take the picture.
        width,height (float): size, in inches, of the output images.
        prefix (str): prefix for the files created. 
     
    Returns: the list of files created (for later removal)
    """
     
    files = []
    ax.figure.set_size_inches(width,height)
     
    for i,angle in enumerate(angles):
     
        ax.view_init(elev = elevation, azim=angle)
        fname = '%s%03d.jpeg'%(prefix,i)
        ax.figure.savefig(fname)
        files.append(fname)
     
    return files

# MP4 file generation [From the web]
def make_movie(files,output, fps=10,bitrate=1800,**kwargs):
    """
    Uses mencoder, produces a .mp4/.ogv/... movie from a list of
    picture files.
    """
     
    output_name, output_ext = os.path.splitext(output)
    command = { '.mp4' : 'mencoder "mf://%s" -mf fps=%d -o %s.mp4 -ovc lavc\
                         -lavcopts vcodec=msmpeg4v2:vbitrate=%d'
                         %(",".join(files),fps,output_name,bitrate)}
                          
    command['.ogv'] = command['.mp4'] + '; ffmpeg -i %s.mp4 -r %d %s'%(output_name,fps,output)
     
    print command[output_ext]
    output_ext = os.path.splitext(output)[1]
    os.system(command[output_ext])
 

# GIF file generation [From the web]
def make_gif(files,output,delay=100, repeat=True,**kwargs):
    """
    Uses imageMagick to produce an animated .gif from a list of
    picture files.
    """
     
    loop = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s'
              %(delay,loop," ".join(files),output))
 
# JPEG strip generation [From the web]
def make_strip(files,output,**kwargs):
    """
    Uses imageMagick to produce a .jpeg strip from a list of
    picture files.
    """
     
    os.system('montage -tile 1x -geometry +0+0 %s %s'%(" ".join(files),output))
     
# Animation [From the web]
def rotanimate(ax, angles, output, **kwargs):
    """
    Produces an animation (.mp4,.ogv,.gif,.jpeg,.png) from a 3D plot on
    a 3D ax
     
    Args:
        ax (3D axis): the ax containing the plot of interest
        angles (list): the list of angles (in degree) under which to
                       show the plot.
        output : name of the output file. The extension determines the
                 kind of animation used.
        **kwargs:
            - width : in inches
            - heigth: in inches
            - framerate : frames per second
            - delay : delay between frames in milliseconds
            - repeat : True or False (.gif only)
    """
         
    output_ext = os.path.splitext(output)[1]
 
    files = make_views(ax,angles, **kwargs)
     
    D = { '.mp4' : make_movie,
          '.ogv' : make_movie,
          '.gif': make_gif ,
          '.jpeg': make_strip,
          '.png':make_strip}
           
    D[output_ext](files,output,**kwargs)
     
    for f in files:
        os.remove(f)
    
# Reversing the color scheme of a given colormap [From the web]
def reverse_colourmap(cmap, name = 'my_cmap_r'):
    """
    In:
    cmap, name
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)
    return my_cmap_r

# Custom GIF producing procedure of the fit model and residual maps
def fits_gif(save_path,output_gif_path,fits_files,cash,psr_alpha,tname,zoom,g1):
 

    gif_model       = []
    gif_resid       = []
       
    smooth    = 0
    
    for fl in range(0,len(fits_files)):

        ra_psr      = 228.4818
        dec_psr     = -59.1358
        glon_psr    = 320.3208
        glat_psr    = -1.1619

        list_model_path       = save_path + fits_files[fl] + '_model.fits'
        list_resid_path       = save_path + fits_files[fl] + '_resid.fits'

        hdu_gif_model    = fits.open(list_model_path)
        hdu_gif_resid    = fits.open(list_resid_path)


        hdu_gam   = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

        header      = hdu_gam[names.hdu_ext].header
        w = wcs.WCS(header)
        w.wcs.ctype = ["RA---GLS", "DEC--GLS"]
        coord1 = ra_psr
        coord2 = dec_psr

        x_psr,y_psr = w.wcs_world2pix(coord1, coord2, names.random_convert_shift)
        
        hdu_gam.close()
        plt.clf()

        reset_display()
        fig = plt.figure()

        ax = fig.add_subplot(111,projection=w)
        ax.patch.set_facecolor('black')

        ax.set_title('Model',fontsize=15)

        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')

        if(zoom.find('all_n') != -1):
            zoom_1    = 160#50
            zoom_2    = 240#350
        if(zoom.find('all_y') != -1):
            zoom_1    = 0
            zoom_2    = 400

        ax.set_xlim(zoom_1,zoom_2)
        ax.set_ylim(zoom_1,zoom_2)        
        
        color_map   = plt.get_cmap('CMRmap')#, 10)
        color_map_2 = plt.get_cmap('seismic')


        #color_map   = reverse_colourmap(color_map)
        #color_map_2 = reverse_colourmap(color_map_2)


        img   = ndimage.gaussian_filter(hdu_gif_model[0].data, sigma=smooth, order=0)


        if(names.gif_norm.find('y') != -1):
            image = ax.imshow(np.sqrt(img), origin="lower", cmap=color_map, vmin=names.cb_extrema[names.analysis,0] , vmax=names.cb_extrema[names.analysis,1])
        else:
            image = ax.imshow(img,origin="lower",cmap=color_map,vmin=img.min(),vmax=img.max())
        #ax.autoscale(False)
        ax.plot(x_psr,y_psr,'+',color='darkgreen',markersize=10,label="PSR B1509-58")

        #plt.colorbar(image,label="Counts")
        plt.colorbar(image, label="Counts ($\sqrt{\;}$ scale)")

        #Convert the parameters for the added component                
        if(tname.find("G") != -1):
            gaus        = []
            gaus.append(g1[fl][0])
            gaus.append(g1[fl][1])
            gaus.append(g1[fl][2])  #/ (2. * np.sqrt(2. * np.log(2.)))
            gaus.append(g1[fl][2] * (1 - g1[fl][3]))
            if(g1[fl][4] >= 0.):
                gaus.append(g1[fl][4] * (180./np.pi))
            if(g1[fl][4] < 0.):
                gaus.append( 360. + (g1[fl][4]* (180./np.pi)))

            gaus.append(g1[fl][5])
            gaus.append(g1[fl][6])


            ellipse = Ellipse(xy=(gaus[0], gaus[1]), width=gaus[2], height=gaus[3], angle = gaus[4], edgecolor='k', fc='None', lw=2)
            ax.add_patch(ellipse)
            ax.plot(gaus[0],gaus[1],'x',color='k',markersize=10)

            x_plane  = np.arange(401)
            x1_plane = 145.74
            x2_plane = 383.74
            y1_plane = 373.41
            y2_plane = 228.31
            y_plane = affine(x_plane,x1_plane,y1_plane,x2_plane,y2_plane)

            if (zoom.find('all_y') != -1):
                ax.plot(x_plane,y_plane,'--',linewidth=3.0,color='c',label="Galactic plane")


        if(tname.find("sh") != -1):
            gaus        = []
            gaus.append(g1[fl][0])            
            gaus.append(g1[fl][1])
            gaus.append(g1[fl][2])
            gaus.append(g1[fl][2] + g1[fl][3])            
            
            gaus.append(g1[fl][4])
            gaus.append(g1[fl][5])
            gaus.append(g1[fl][6])

            circle_in  = Circle((gaus[0],gaus[1]),radius=gaus[2],fc='None',edgecolor='k',lw=1,alpha=0.8)
            circle_out = Circle((gaus[0],gaus[1]),radius=gaus[3],fc='None',edgecolor='k',lw=1)            

            ax.add_patch(circle_in)
            ax.add_patch(circle_out)
            ax.plot(gaus[0],gaus[1],'x',color='k',markersize=10)             
            

        if(tname.find("dc") != -1):  
            gaus       = []
            gaus.append(g1[fl][0])
            gaus.append(g1[fl][1])
            gaus.append(g1[fl][2])
            gaus.append(g1[fl][3])
            gaus.append(g1[fl][4])
            gaus.append(g1[fl][5])
            gaus.append(g1[fl][6])

            circle = Circle((gaus[0],gaus[1]),radius=gaus[2],fc='k',edgecolor='k',lw=3,alpha=0.3)
            
            ax.add_patch(circle)
            ax.plot(gaus[0],gaus[1],'x',color='k',ms=10)

            
        if(tname.find("dc") == -1 and tname.find("sh") == -1 and tname.find("G") == -1):  
            gaus        = np.zeros(6)
            gaus[5]     = (g1[fl][0])
            

        if(tname.find("dc") != -1 or tname.find("sh") != -1 or tname.find("G") != -1):

            if(cash[fl] == min(cash)):

                if(tname.find("G") != -1):

                    binsize = 0.01

                    ax.annotate(r"$\alpha$       : " + str(psr_alpha[fl]), (zoom_1 + 5, zoom_1 + 5),fontsize=20, color='w')
                    #ax.annotate(r"$\alpha$       : " + str(psr_alpha[fl]) + "\nCASH : " + str(cash[fl]) + "\nGaussian A / X-ray A :" + str((gaus[6] * 2. * np.pi * (gaus[2] / (2. * np.sqrt(2. * np.log(2.)))) ** 2) / gaus[5])[:5],(zoom_1 + 5, zoom_1 + 5), fontsize=15, color='w')
                    #ax.annotate("a         : " + str(psr_alpha[fl]) + "\nCASH : " + str(cash[fl]) + "\nX-ray A :" + str(gaus[5]) + "\nGaussian A / X-ray A :" + str( (gaus[6] * 2. * np.pi * fwhm_sigma(gaus[2],binsize)**2) / gaus[5])[:5], (zoom_1 + 5, zoom_1 + 5),fontsize=15, color='r')
                else :
                    ax.annotate(r"$\alpha$       : " + str(psr_alpha[fl]), (zoom_1 + 5, zoom_1 + 5), fontsize=20, color='w')
                    #ax.annotate("a         : " + str(psr_alpha[fl]) + "\nCASH : " + str(cash[fl]) + "\nX-ray A :" + str(gaus[5]) + "\nAmpl / X-ray A :" + str( (gaus[6] ) / gaus[5])[:7], (zoom_1 + 5, zoom_1 + 5),fontsize=15, color='r')

        else:
            ax.annotate(r"$\alpha$       : " + str(psr_alpha[fl]),(zoom_1 + 5, zoom_1 + 5), fontsize=20, color='w')
            #ax.annotate("a         : "+str(psr_alpha[fl])+"\nCASH : "+str(cash[fl])+"\nX-ray A :"+str(gaus[5]),(zoom_1+5,zoom_1+5),fontsize=15,color='w')
            if(cash[fl] == min(cash)):
                ax.annotate(r"$\alpha$       : " + str(psr_alpha[fl]),(zoom_1 + 5, zoom_1 + 5), fontsize=20, color='w')
                #ax.annotate("a         : "+str(psr_alpha[fl])+"\nCASH : "+str(cash[fl])+"\nX-ray A :"+str(gaus[5]),(zoom_1+5,zoom_1+5),fontsize=15,color='r')
                    
            
        ax.tick_params(axis='both',which='both',labelsize=20,color='w') 
        ax.legend(fontsize=15,numpoints=1,loc=2)
        plt.savefig(output_gif_path+fits_files[fl]+'_model_'+zoom+'.png',dpi=100)        
        plt.clf()

        reset_display()
        fig = plt.figure()

        ax = fig.add_subplot(111,projection=w)
        ax.patch.set_facecolor('black')


        ax.set_title('Residuals',fontsize=15)
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')

        if (zoom.find('all_n') != -1):
            zoom_1 = 160
            zoom_2 = 240
        if (zoom.find('all_y') != -1):
            zoom_1 = 0
            zoom_2 = 400

        ax.set_xlim(zoom_1,zoom_2)
        ax.set_ylim(zoom_1,zoom_2)        
        
        img = ndimage.gaussian_filter(hdu_gif_resid[0].data, sigma=smooth+1, order=0)

        if (names.gif_norm.find('y') != -1) :
            image = ax.imshow(img, origin="lower", cmap=color_map_2, vmin=names.cb_resids[names.analysis,0], vmax=names.cb_resids[names.analysis,1])
        else:
            image = ax.imshow(img, origin="lower", cmap=color_map, vmin=img.min(), vmax=img.max())

        ax.plot(x_psr,y_psr,'+',color='darkgreen',markersize=10,label="PSR B1509-58")        

        plt.colorbar(image,label="Counts")        

        # ax.annotate("a         : "+str(psr_alpha[fl])+"\nCASH : "+str(cash[fl]),(zoom_1+5,zoom_1+5),fontsize=15,color='w')
        # if(cash[fl] == min(cash)):
        #     ax.annotate("a         : "+str(psr_alpha[fl])+"\nCASH : "+str(cash[fl]),(zoom_1+5,zoom_1+5),fontsize=15,color='b')

        ax.tick_params(axis='both',which='both',labelsize=20,color='w') 
        ax.legend(fontsize=15,numpoints=1,loc=1)
        plt.savefig(output_gif_path + fits_files[fl]+'_resid_'+zoom+'.png',dpi=100)        
        plt.clf()

        
        hdu_gif_model.close()
        hdu_gif_resid.close()
        
        gif_model.append(output_gif_path + fits_files[fl]+'_model_'+zoom+'.png')
        gif_resid.append(output_gif_path + fits_files[fl]+'_resid_'+zoom+'.png')
                
    make_gif(files=gif_model,delay=100,output=save_path + tname+"_model_"+zoom+".gif")
    make_gif(files=gif_resid,delay=100,output=save_path + tname+"_resid_"+zoom+".gif")

# Custom PNG producing procedure of the fit model and residual maps
def fits_png(save_path, output_gif_path, fits_files, cash, psr_alpha, tname, zoom, g1, cb_reverse):

    smooth = 0

    ra_psr = 228.4818
    dec_psr = -59.1358
    glon_psr = 320.3208
    glat_psr = -1.1619

    list_model_path = save_path + fits_files + '_model.fits'
    list_resid_path = save_path + fits_files + '_resid.fits'

    hdu_gif_model = fits.open(list_model_path)
    hdu_gif_resid = fits.open(list_resid_path)

    hdu_gam = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

    header = hdu_gam[names.hdu_ext].header
    w = wcs.WCS(header)
    w.wcs.ctype = ["RA---GLS", "DEC--GLS"]
    coord1 = ra_psr
    coord2 = dec_psr

    x_psr, y_psr = w.wcs_world2pix(coord1, coord2, names.random_convert_shift)

    hdu_gam.close()
    plt.clf()

    reset_display()
    fig = plt.figure()

    ax = fig.add_subplot(111, projection=w)
    ax.patch.set_facecolor('black')

    ax.set_title('Model', fontsize=15)

    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')

    if (zoom.find('all_n') != -1):
        zoom_1 = 160  # 50
        zoom_2 = 240  # 350
    if (zoom.find('all_y') != -1):
        zoom_1 = 0
        zoom_2 = 400

    ax.set_xlim(zoom_1, zoom_2)
    ax.set_ylim(zoom_1, zoom_2)


    color_map = plt.get_cmap('CMRmap')  # , 10)
    color_map_2 = plt.get_cmap('seismic')

    if(cb_reverse.find('y') != -1):
        color_map   = reverse_colourmap(color_map)
        color_map_2 = reverse_colourmap(color_map_2)


    img = ndimage.gaussian_filter(hdu_gif_model[0].data, sigma=smooth, order=0)

    if (names.gif_norm.find('y') != -1):
        image = ax.imshow(np.sqrt(img), origin="lower", cmap=color_map, vmin=names.cb_extrema[names.analysis, 0],
                              vmax=names.cb_extrema[names.analysis, 1])
    else:
        image = ax.imshow(img, origin="lower", cmap=color_map, vmin=img.min(), vmax=img.max())
        # ax.autoscale(False)
    ax.plot(x_psr, y_psr, '+', color='darkgreen', markersize=10, label="PSR B1509-58")

    plt.colorbar(image, label="Counts ($\sqrt{\;}$ scale)")

    # Convert the parameters for the added component
    if (tname.find("G") != -1):
        gaus = []
        gaus.append(g1[0])
        gaus.append(g1[1])
        gaus.append(g1[2])  # / (2. * np.sqrt(2. * np.log(2.)))
        gaus.append(g1[2] * (1 - g1[3]))
        if (g1[4] >= 0.):
            gaus.append(g1[4] * (180. / np.pi))
        if (g1[4] < 0.):
            gaus.append(360. + (g1[4] * (180. / np.pi)))

        gaus.append(g1[5])
        gaus.append(g1[6])

        ellipse = Ellipse(xy=(gaus[0], gaus[1]), width=gaus[2], height=gaus[3], angle=gaus[4], edgecolor='w',
                              fc='None', lw=2)
        ax.add_patch(ellipse)
        ax.plot(gaus[0], gaus[1], 'x', color='k', markersize=10)

        if (cb_reverse.find('y') != -1):
            ellipse = Ellipse(xy=(gaus[0], gaus[1]), width=gaus[2], height=gaus[3], angle=gaus[4], edgecolor='k',
                              fc='None', lw=2)
            ax.add_patch(ellipse)
            ax.plot(gaus[0], gaus[1], 'x', color='w', markersize=10)

    x_plane = np.arange(401)
    x1_plane = 145.74
    x2_plane = 383.74
    y1_plane = 373.41
    y2_plane = 228.31
    y_plane = affine(x_plane, x1_plane, y1_plane, x2_plane, y2_plane)

    if (zoom.find('all_y') != -1):
        ax.plot(x_plane, y_plane, '--', linewidth=3.0, color='c', label="Galactic plane")

    if (tname.find("sh") != -1):
        gaus = []
        gaus.append(g1[0])
        gaus.append(g1[1])
        gaus.append(g1[2])
        gaus.append(g1[2] + g1[3])

        gaus.append(g1[4])
        gaus.append(g1[5])
        gaus.append(g1[6])

        circle_in = Circle((gaus[0], gaus[1]), radius=gaus[2], fc='None', edgecolor='w', lw=1, alpha=0.8)
        circle_out = Circle((gaus[0], gaus[1]), radius=gaus[3], fc='None', edgecolor='w', lw=1)
        ax.add_patch(circle_in)
        ax.add_patch(circle_out)
        ax.plot(gaus[0], gaus[1], 'x', color='k', markersize=10)

        if (cb_reverse.find('y') != -1):
            circle_in = Circle((gaus[0], gaus[1]), radius=gaus[2], fc='None', edgecolor='k', lw=1, alpha=0.8)
            circle_out = Circle((gaus[0], gaus[1]), radius=gaus[3], fc='None', edgecolor='k', lw=1)
            ax.add_patch(circle_in)
            ax.add_patch(circle_out)
            ax.plot(gaus[0], gaus[1], 'x', color='k', markersize=10)

    if (tname.find("dc") != -1):
        gaus = []
        gaus.append(g1[0])
        gaus.append(g1[1])
        gaus.append(g1[2])
        gaus.append(g1[3])
        gaus.append(g1[4])
        gaus.append(g1[5])
        gaus.append(g1[6])

        circle = Circle((gaus[0], gaus[1]), radius=gaus[2], fc='k', edgecolor='w', lw=3, alpha=0.3)
        ax.add_patch(circle)
        ax.plot(gaus[0], gaus[1], 'x', color='k', ms=10)

        if (cb_reverse.find('y') != -1):
            circle = Circle((gaus[0], gaus[1]), radius=gaus[2], fc='k', edgecolor='k', lw=3, alpha=0.3)
            ax.add_patch(circle)
            ax.plot(gaus[0], gaus[1], 'x', color='k', ms=10)

    if (tname.find("dc") == -1 and tname.find("sh") == -1 and tname.find("G") == -1):
        gaus = np.zeros(6)
        gaus[5] = (g1[0])

    ax.annotate(r"$\alpha$       : " + str(psr_alpha), (zoom_1 + 5, zoom_1 + 5), fontsize=20, color='w')

    if(cb_reverse.find('y') != -1):
        ax.annotate(r"$\alpha$       : " + str(psr_alpha), (zoom_1 + 5, zoom_1 + 5), fontsize=20, color='k')



    ax.tick_params(axis='both', which='both', labelsize=20, color='w')
    ax.legend(fontsize=15, numpoints=1, loc=2)

    plt.savefig(output_gif_path + fits_files[:18] + tname + '_model_' + zoom + cb_reverse + '.png', dpi=100)
    plt.clf()

    reset_display()
    fig = plt.figure()

    ax = fig.add_subplot(111, projection=w)
    ax.patch.set_facecolor('black')

    ax.set_title('Residuals', fontsize=15)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')

    if (zoom.find('all_n') != -1):
        zoom_1 = 160
        zoom_2 = 240
    if (zoom.find('all_y') != -1):
        zoom_1 = 0
        zoom_2 = 400

    ax.set_xlim(zoom_1, zoom_2)
    ax.set_ylim(zoom_1, zoom_2)

    img = ndimage.gaussian_filter(hdu_gif_resid[0].data, sigma=smooth + 1, order=0)

    if (names.gif_norm.find('y') != -1):
        image = ax.imshow(img, origin="lower", cmap=color_map_2, vmin=names.cb_resids[names.analysis, 0],
                              vmax=names.cb_resids[names.analysis, 1])
    else:
        image = ax.imshow(img, origin="lower", cmap=color_map, vmin=img.min(), vmax=img.max())

    ax.plot(x_psr, y_psr, '+', color='darkgreen', markersize=10, label="PSR B1509-58")

    plt.colorbar(image, label="Counts")


    ax.tick_params(axis='both', which='both', labelsize=20, color='w')
    ax.legend(fontsize=15, numpoints=1, loc=1)
    plt.savefig(output_gif_path + fits_files[:18] + tname + '_resid_' + zoom + cb_reverse + '.png', dpi=100)
    plt.clf()

    hdu_gif_model.close()
    hdu_gif_resid.close()

# Plotting the energy slicing for a _ vs energy plot/histogram
def slice_energy_plot(E_cuts, colors):

    bins   = len(E_cuts) - 1

    for i in range(0, bins):
        plt.axvspan(E_cuts[i], E_cuts[i+1], color=colors[i], linestyle='solid', alpha=0.2)

    return

# Plotting the excess label for each energy bin for an excess vs energy plot/histogram
def slice_label_excess_plot(E_cuts, Ex_slices, X_gamma, colors):

    for i in range(0,len(Ex_slices) -1):

        diff_exc = Ex_slices[i+1] - Ex_slices[i]

        plt.annotate(str(format(diff_exc, '.1f')), ((E_cuts[i] + E_cuts[i+1]) / 2.3, max(X_gamma) + 3.), fontsize=12, color=colors[i])

    return