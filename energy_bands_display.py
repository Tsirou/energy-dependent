import os
import names

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import scipy.ndimage as ndimage


from astropy import wcs
from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.patches import Ellipse,Circle

from conversion import affine,progress,radial_profile,crop,crop_circle,convert_image
from display import reset_display,make_gif



fits_files       = []
zoom             = "all_n"

#labels           = [" > 0.3 TeV", "[0.3 - 0.6] TeV", "[0.6 - 0.9] TeV", "[0.9 - 3.0] TeV", " > 3.0 TeV"]
labels           = ["[0.3 - 0.6] TeV", "[0.6 - 1.2] TeV", "[1.2 - 2.4] TeV", " > 2.4 TeV"]

save_path        = names.path_template[0]
output_gif_path  = names.path[0]

# for f in range (0,len(names.filename_gamma)):
#     fits_files.append(names.filename_gamma[f])

for f in range (0,len(names.filename_excess) - 1):
    fits_files.append(names.filename_excess[f])



g1               = []



def fits_gif(save_path, output_gif_path, fits_files, zoom, g1):
    gif_model = []

    smooth = 1

    for fl in range(1, len(fits_files)):

        progress(fl,len(fits_files),"\n"+labels[fl])

        ra_psr   = 228.4818
        dec_psr  = -59.1358
        glon_psr = 320.3208
        glat_psr = -1.1619

        list_model_path = save_path + fits_files[fl]

        hdu_gif_model   = fits.open(list_model_path)

        hdu_gam         = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

        header          = hdu_gam[names.hdu_ext].header
        w               = wcs.WCS(header)
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

        ax.set_title('MSH 15-5$\mathit{2}$ energy-dependent analyses', fontsize=13)

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

        # color_map   = reverse_colourmap(color_map)
        # color_map_2 = reverse_colourmap(color_map_2)


        img = ndimage.gaussian_filter(hdu_gif_model[names.hdu_ext].data, sigma=smooth, order=   0)

        # if (names.gif_norm.find('y') != -1):
        #     image = ax.imshow(img, origin="lower", cmap=color_map, vmin=1.0, vmax=4.0)
        # else:
        image = ax.imshow(img, origin="lower", cmap=color_map, vmin=img.min(), vmax=img.max())

        ax.plot(x_psr, y_psr, '+', color='darkgreen', markersize=12, label="PSR B1509-58")

        plt.colorbar(image, label="Counts")

        # # Convert the parameters for the added component
        # gaus = []
        # gaus.append(g1[fl][0])
        # gaus.append(g1[fl][1])
        # gaus.append(g1[fl][2])  # / (2. * np.sqrt(2. * np.log(2.)))
        # gaus.append(g1[fl][2] * (1 - g1[fl][3]))
        # if (g1[fl][4] >= 0.):
        #     gaus.append(g1[fl][4] * (180. / np.pi))
        # if (g1[fl][4] < 0.):
        #     gaus.append(360. + (g1[fl][4] * (180. / np.pi)))
        #
        # gaus.append(g1[fl][5])
        # gaus.append(g1[fl][6])
        #
        # ellipse = Ellipse(xy=(gaus[0], gaus[1]), width=gaus[2], height=gaus[3], angle=gaus[4], edgecolor='k',
        #                       fc='None', lw=2)
        # ax.add_patch(ellipse)
        # ax.plot(gaus[0], gaus[1], 'x', color='k', markersize=10)

        x_plane = np.arange(401)
        x1_plane = 145.74
        x2_plane = 383.74
        y1_plane = 373.41
        y2_plane = 228.31
        y_plane = affine(x_plane, x1_plane, y1_plane, x2_plane, y2_plane)

        if (zoom.find('all_y') != -1):
            ax.plot(x_plane, y_plane, '--', linewidth=3.0, color='c', label="Galactic plane")

            ax.annotate(labels[fl], (zoom_1 + 10, zoom_2 - 84), fontsize=18, color='w')
            ax.legend(fontsize=15, numpoints=1, loc=4)


        else:
            ax.annotate(labels[fl], (zoom_1 + 5, zoom_1 +5), fontsize=18, color='w')
            ax.legend(fontsize=15, numpoints=1, loc=1)



        ax.tick_params(axis='both', which='both', labelsize=20, color='w')
        plt.savefig(output_gif_path + fits_files[fl] + '_energy-dependent_' + zoom + '.png', dpi=100)
        plt.clf()


        hdu_gif_model.close()

        gif_model.append(output_gif_path + fits_files[fl] + '_energy-dependent_' + zoom + '.png')

    make_gif(files=list(reversed(gif_model)), delay=120, output=output_gif_path + "energy-dependent_excess_r_" + zoom + ".gif")
    #make_gif(files=gif_model, delay=120, output=output_gif_path + "energy-dependent_excess_" + zoom + ".gif")

    return


def profile(save_path, fits_files, output_gif_path):
    ra_psr = 228.4818
    dec_psr = -59.1358
    glon_psr = 320.3208
    glat_psr = -1.1619

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

    plt.style.use('dark_background')

    fig = plt.figure()

    fig.set_size_inches(10.5,8.5)

    ax = fig.add_subplot(111)
    ax.patch.set_facecolor('black')
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')

    ax.xaxis.label.set_color('w')

    ax.set_title('MSH 15-5$\mathit{2}$ energy-dependent profiles', fontsize=13, color='w')

    ax.set_xlabel('Distance to PSR B1509-58')
    ax.set_ylabel('Excess')

    colors = ["yellow","cornflowerblue", "springgreen", "orange", "indianred"]

    for fl in range(1,len(fits_files)):

        list_model_path = save_path + fits_files[fl]

        hdu_gif_model = fits.open(list_model_path)

        source        = hdu_gif_model[names.hdu_ext].data

        img = ndimage.gaussian_filter(source, sigma=0, order=0)

        cut  = radial_profile(source, [x_psr,y_psr])


        x_plane  = np.arange(401)
        x1_plane = 145.74
        x2_plane = 383.74
        y1_plane = 373.41
        y2_plane = 228.31
        y_plane = affine(x_plane, x1_plane, y1_plane, x2_plane, y2_plane)

        x_cut   = np.arange(401)
        y_cut   = affine(x_cut, x_psr, 0., y_psr, 0.)


        cuty = y_cut[:401]
        cuts = []
        for i in range(0,len(source)):
            for j in range(0,len(source)):
                if( abs( (cuty[i] - j)) < 2):
                    cuts.append(source[i,j])



        ax.plot(cut[:40],color=colors[fl],label=labels[fl],marker="+",markersize=12, linestyle='--')
        #ax.plot(x_cut,cut, color=colors[fl], label=labels[fl], marker="+", markersize=12, linestyle='--')

    ax.legend(fontsize=15, numpoints=1, loc=1)

    ax.tick_params(axis='both', which='both', labelsize=20, color='w')

    plt.show()

    plt.savefig(output_gif_path + '_energy-dependent_profiles' + zoom + '.png', dpi=100)
    plt.clf()

    hdu_gif_model.close()


    return


def fits_cuts(save_path, output_gif_path, fits_files,zoom):
    gif_model = []

    smooth = 1

    #colors = ["yellow","cornflowerblue", "springgreen", "orange", "indianred"]
    colors = ["cornflowerblue", "springgreen", "orange", "indianred"]

    for fl in range(0, len(fits_files)):

        progress(fl, len(fits_files), "\n" + labels[fl])

        ra_psr  = 228.4818
        dec_psr = -59.1358

        list_model_path = save_path + fits_files[fl]

        hdu_gif_model = fits.open(list_model_path)
        source        = hdu_gif_model[names.hdu_ext].data
        hdu_gam       = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

        header = hdu_gam[names.hdu_ext].header
        ww = wcs.WCS(header)
        ww.wcs.ctype = ["RA---GLS", "DEC--GLS"]
        coord1 = ra_psr
        coord2 = dec_psr

        x_psr, y_psr = ww.wcs_world2pix(coord1, coord2, names.random_convert_shift)

        x_pix_psr = int(x_psr)
        y_pix_psr = int(y_psr)

        hdu_gam.close()
        plt.clf()

        reset_display()
        fig = plt.figure()
        fig.set_size_inches(12.5, 10.5)
        plt.style.use('dark_background')

        if (zoom.find('all_n') != -1):
            zoom_1 = 160
            zoom_2 = 240

        source_i = source
        source = crop(source,zoom_2-zoom_1)

        h, w = source.shape
        gs = gridspec.GridSpec(2, 2, width_ratios=[w, w * .2], height_ratios=[h, h * .2])

        ax = [plt.subplot(gs[0], projection=ww), ]
        ax.append(plt.subplot(gs[1], sharey=ax[0]))
        ax.append(plt.subplot(gs[2], sharex=ax[0]))
        ax.append(plt.subplot(gs[3]))


        ax[0].set_xlim(zoom_1, zoom_2)
        ax[0].set_ylim(zoom_1, zoom_2)

        ax[0].coords.grid(True, color='k', ls='solid')

        ax[0].coords[0].set_axislabel('RA')
        ax[0].coords[1].set_axislabel('Dec')

        overlay = ax[0].get_coords_overlay('galactic')
        overlay.grid(color='grey', ls='dotted')
        overlay[0].set_axislabel('l')
        overlay[1].set_axislabel('b')

        color_map = plt.get_cmap('seismic')

        img = ndimage.gaussian_filter(source_i, sigma=smooth, order=0)
        image = ax[0].imshow(img, origin="lower", cmap=color_map, vmin=img.min(), vmax=img.max())
        #image = ax[0].imshow(img, origin="lower", cmap=color_map, vmin=0., vmax=1.2)

        cax   = fig.add_axes([0.125, 0.93, 0.585, 0.03 ])
        cbar = fig.colorbar(image, cax=cax,orientation='horizontal')
        cbar.set_label(label='Excess', size=15, weight='bold', horizontalalignment='right')

        ax[0].plot(x_psr, y_psr, '+', color='yellow', markersize=15, label="PSR B1509-58", linewidth=10)
        ax[0].axhline(y=y_psr, xmin=0, xmax=len(source), c="yellow", linewidth=2.0, linestyle='--')
        ax[0].axvline(x=x_psr, ymin=0, ymax=len(source), c="yellow", linewidth=2.0, linestyle='--')

        sliced_x = np.arange(400)
        sliced_y = affine(sliced_x,x_psr,y_psr,zoom_1,zoom_1)
        ax[0].plot(sliced_x[zoom_1 + 15:zoom_2 -15],sliced_y[zoom_1 + 15:zoom_2 - 15],'--',color='greenyellow',linewidth=3.0)

        rect = patches.Rectangle((zoom_1 + 2, zoom_1 + 2),30,7,facecolor='k',alpha=1.0)
        box  = patches.Rectangle((zoom_1 + 1, zoom_1 + 1),32,9,facecolor='w',alpha=1.0)

        ax[0].annotate(labels[fl], (zoom_1 + 5, zoom_1 + 5), fontsize=18, color=colors[fl])
        ax[0].add_patch(box)
        ax[0].add_patch(rect)
        ax[0].legend(fontsize=15, numpoints=1, loc=1)
        ax[0].tick_params(axis='both', which='both', labelsize=20, color='w')

        x, y = np.linspace(zoom_1, zoom_2, zoom_2-zoom_1), np.linspace(zoom_1, zoom_2, zoom_2-zoom_1)
        X, Y = np.meshgrid(x, y)

        #center = [w/2, h/2]
        center = [int(x_pix_psr * ( (zoom_2-zoom_1) / 400.)),int(y_pix_psr * ( (zoom_2-zoom_1) / 400.))]

        ax[1].plot(source[:, center[0]], Y[:, (zoom_2-zoom_1) / 2], '+' , source[:,  center[0]], Y[:, (zoom_2-zoom_1) / 2], color=colors[fl])
        ax[1].axis([source[:, center[0]].max(), source[:, center[0]].min(), Y.min(), Y.max()])
        ax[1].axhline(y=y_psr, xmin=0, xmax=len(source), c="yellow", linewidth=2.0, linestyle='--')
        ax[1].yaxis.set_label_position("right")
        ax[1].yaxis.tick_right()
        ticks = ax[1].get_yticks()
        xpsr  = [x_psr for xp in ticks]
        cel1, cel2 = ww.wcs_pix2world(xpsr, ticks, names.random_convert_shift)


        ticks_cel = []
        for e in range(0,len(cel2)):
            temp_val   = (dec_psr - cel2[e]) * 60.
            ticks_cel.append( round(temp_val,2))

        ax[1].set_yticklabels(ticks_cel)
        ax[1].set_ylabel('Distance to the pulsar (arcmin)')
        ax[1].set_xlabel('Excess')

        ax[2].plot(X[(zoom_2-zoom_1) / 2, :], source[center[1], :], '+', X[(zoom_2-zoom_1) / 2, :], source[center[1], :], color=colors[fl])
        ax[2].axvline(x=x_psr, ymin=0, ymax=len(source), c="yellow", linewidth=2.0, linestyle='--')

        ticks = ax[2].get_xticks()
        ypsr  = [y_psr for xp in ticks]
        cel1, cel2 = ww.wcs_pix2world(ticks, ypsr, names.random_convert_shift)


        ticks_cel = []
        for e in range(0,len(cel1)):
            temp_val   = (ra_psr - cel1[e]) * 60.
            ticks_cel.append(round(temp_val,2))

        fig.canvas.draw()
        labels_t = [item.get_text() for item in ax[2].get_xticklabels()]
        labels_x = []
        labels_x.append(u'')
        for ele in ticks_cel[::2]:
            labels_x.append(str(ele))
        labels_x.append(u'')

        labels_t = labels_x
        ax[2].set_xticklabels(labels_t)
        ax[2].set_xlabel('Distance to the pulsar (arcmin)')
        ax[2].set_ylabel('Excess')

        cut_prof = radial_profile(source_i, [x_psr, y_psr])
        ax[3].plot(cut_prof[:50],'+',color=colors[fl])
        ax[3].set_ylim(0.,2.0)
        ax[3].set_xlim(0.,40.)
        ax[3].axvline(x=0., ymin=0, ymax=4.0, c="greenyellow", linewidth=2.0, linestyle='--')

        ticks = ax[3].get_xticks()

        ticks = ticks * 0.01 * 60
        ax[3].set_xticklabels(ticks)
        ax[3].set_xlabel('Distance to the pulsar (arcmin)')
        ax[3].set_ylabel('Excess')


        plt.savefig(output_gif_path + fits_files[fl] + '_energy_distribution_2d0_ebins_hires' + zoom + '.png', dpi=100)
        plt.clf()

        hdu_gif_model.close()

        gif_model.append(output_gif_path + fits_files[fl] + '_energy_distribution_2d0_ebins_hires' + zoom + '.png')

    make_gif(files=list(reversed(gif_model)), delay=160,
             output=output_gif_path + "energy-dependent_excess_distribution_2d0_ebins_hires" + zoom + ".gif")


    return



################################################################################################
################################################################################################

#fits_gif(save_path, output_gif_path,fits_files,zoom,g1)
#profile(save_path,fits_files,output_gif_path)
fits_cuts(save_path,output_gif_path,fits_files,zoom)
