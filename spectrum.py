
#Stupid old stuff
# import os
# import names
#
# from PyAstronomy import pyasl
# from astropy.io import fits
#
#
#
# # c          = 2.9992458e8   # m.s-1
# # E_nu       = 1.0e12        # eV
# #
# # wavelength = h * c /E_nu / 1.0e-10
# #
# # print wavelength
# #
# # os._exit(0)
# #
# # hdu_g   = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])
# #
# # image   = hdu_g[names.hdu_ext].data
# #
# # for i in range(0,len(image)):
# #     for j in range(0,len(image)):
# #         photons = image[i,j]
# #         flux    = pyasl.photons2flux(wavelength,photons)
# #
# # print(flux)
# #
# # hdu_g.close()
#
#
#
#
# for comp in range(0,3):
#     set_full_model(comp, bkg + psf((gauss2d.g0 + gauss2d.g1 + gauss2d.g2)) * exposure)
#
#     if(comp == 0 ):
#         g0.integrate=True
#         g1.integrate=False
#         g2.integrate=False
#
#     if(comp == 1 ):
#         g0.integrate=True
#         g1.integrate=True
#         g2.integrate=False
#
#     if(comp == 2 ):
#         g0.integrate=True
#         g1.integrate=True
#         g2.integrate=True

#################################################################

import names

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from astropy.coordinates import SkyCoord, Angle
from gammapy import datasets
from gammapy.data import DataStore
from gammapy.image import SkyImage, SkyMask, SkyImageList
from gammapy.scripts import StackedObsImageMaker, SingleObsImageMaker, SpectrumObservation
from gammapy.utils.energy import Energy, EnergyBounds
from gammapy.catalog import SourceCatalogGammaCat, SourceCatalogHGPS
from gammapy.cube import exposure_cube, SkyCube, StackedObsCubeMaker
from gammapy.irf import EffectiveAreaTable2D
from gammapy.spectrum import models

from regions import CircleSkyRegion

from conversion import progress, pixarea_correction
from astropy.io import fits

######################################################################################
# Script for generating exposure maps with gammapy for the energy dependent analyses #
######################################################################################

#This code works if I twist a bit the error messages when reading the header of the fits files. Have to be careful!

num_slices         = 10

path               = names.path_g[names.analysis]

gammacat           = SourceCatalogGammaCat()
gammacat_source    = gammacat['MSH 15-52']

dataset_directory  = names.path_data + "release-1.0/Model_Deconvoluted_Prod26/Mpp_Std/"
data_store         = DataStore.from_dir(dataset_directory)

plt.style.use('ggplot')


obs_ids            = np.loadtxt(names.path_runs[names.analysis]+names.filename_runs[names.analysis],dtype=int,skiprows=1,unpack=True)

source_pos          = SkyCoord(228.5290, -59.1575,unit='deg')

on_region = CircleSkyRegion(source_pos, 0.3 * u.deg)


ref_image = SkyImage.empty(
    nxpix=400, nypix=400, binsz=0.01,
    xref=source_pos.ra.deg, yref=source_pos.dec.deg,
    coordsys='CEL', proj='GLS'
)


energy_band  = Energy([names.analysis_energy_low[names.analysis], names.analysis_energy_high[names.analysis]], 'TeV')
energy_range = [Energy(names.analysis_energy_low[names.analysis], "TeV"), Energy(names.analysis_energy_high[names.analysis], "TeV"), num_slices]

offset_band  = Angle([0, 3.0], 'deg')
offset_range = [Angle(0.0, "deg"), Angle(3.0 , "deg"), num_slices]


exclusion_mask = SkyMask.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
# exclusion_mask = exclusion_mask.reproject(reference=ref_image)

counts_image       = SkyImage.empty_like(ref_image)
for obs in range(0,len(obs_ids)-3):

    events = data_store.obs(obs_id=obs_ids[obs]).events
    counts_image.fill_events(events)

hdu_index    = ['events', 'gti', 'aeff_2d', 'edisp_2d', 'psf_table', 'psf_3gauss', 'psf_king', 'bkg_2d', 'bkg_3d']
hdu_index_s  = ['events', 'aeff', 'edisp', 'psf_king']

events_list  = data_store.obs_list(obs_id=obs_ids[:-3])
events_data  = data_store.copy_obs(obs_id=obs_ids[:-3], outdir=path+'selection/', clobber=True)

events_store = DataStore.from_dir(path + '/selection/')


images            = SkyImageList()
spectral_index    = 2.2
for_integral_flux = False


model = models.ExponentialCutoffPowerLaw(
    index=1.9,
    lambda_=0.16,
    amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference = 1 * u.TeV,
)

flux_point_binning = EnergyBounds.equal_log_spacing(0.7, 30, 5, u.TeV)


config = dict(
    outdir = None,
    background = dict(
        on_region=on_region,
        exclusion_mask=exclusion_mask,
        min_distance = 0.1 * u.rad,
    ),
    extraction = dict(containment_correction=False),
    fit = dict(
        model=model,
        stat='wstat',
        forward_folded=True,
        fit_range = flux_point_binning[[0, -1]]
    ),
    fp_binning=flux_point_binning
)



ana = SpectrumObservation(
    observations=events_list,
    config=config,
)
ana.run()




