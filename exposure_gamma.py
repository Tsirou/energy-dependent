import names

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from astropy.coordinates import SkyCoord, Angle
from gammapy import datasets
from gammapy.data import DataStore
from gammapy.image import SkyImage, SkyMask, SkyImageList
from gammapy.scripts import StackedObsImageMaker, SingleObsImageMaker
from gammapy.utils.energy import Energy
from gammapy.catalog import SourceCatalogGammaCat, SourceCatalogHGPS
from gammapy.cube import exposure_cube, SkyCube, StackedObsCubeMaker
from gammapy.irf import EffectiveAreaTable2D

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


ref_image = SkyImage.empty(
    nxpix=400, nypix=400, binsz=0.01,
    xref=source_pos.ra.deg, yref=source_pos.dec.deg,
    coordsys='CEL', proj='GLS'
)


energy_band  = Energy([names.analysis_energy_low[names.analysis], names.analysis_energy_high[names.analysis]], 'TeV')
energy_range = [Energy(names.analysis_energy_low[names.analysis], "TeV"), Energy(names.analysis_energy_high[names.analysis], "TeV"), num_slices]

offset_band  = Angle([0, 3.0], 'deg')
offset_range = [Angle(0.0, "deg"), Angle(3.0 , "deg"), num_slices]


# exclusion_mask = SkyMask.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
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

total_exposure = SkyImage.empty_like(ref_image, name='exposure')

observation    = 0
for obs_id in events_store.obs_table['OBS_ID']:

    observation = observation + 1
    progress(observation, len(obs_ids[:-3])," loading ")

    obs = events_store.obs(obs_id)
    obs_image = SingleObsImageMaker(obs, ref_image, energy_band, offset_band)

    if len(obs_image.events.table) <= 0:
        continue
    else:
        obs_image.exposure_image(spectral_index, for_integral_flux)
        total_exposure.data += obs_image.images['exposure'].data


images['exposure'] = total_exposure


plt.clf()
images['exposure'].plot(add_cbar=True,cmap='viridis')
images['exposure'].write(names.path_expo[names.analysis] + 'exposure_map_'+names.tags[names.analysis]+'.fits', overwrite=True)

plt.clf()


# PROSOXI : should change names.filename_expo for the exposure maps without the applied correction in the names.py file!!!

hduexp             = fits.open(names.path_expo[names.analysis] + names.filename_expo[names.analysis])

dimension = len(hduexp[0].data)

corr_exposure_fits = np.ndarray((dimension, dimension))

corr_exposure_fits = pixarea_correction(hduexp[0].data, names.ra_pt, names.dec_pt, names.res_ana, corr_exposure_fits)

hduexp[0].data   = corr_exposure_fits[:,:]

hduexp.writeto(names.path_expo[names.analysis] + 'exposure_map_' + names.tags[names.analysis] + '_PApix' + '.fits.gz', overwrite=True)

hduexp.close()

