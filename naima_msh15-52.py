import naima
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

from astropy.io import ascii
from astropy.table import Table


path   = "/home/tsirou/Documents/Analyses/MSH_15-52/energy-dependent/MultiEnergyBins/spectra/"

#data_pa_ECPL  = ascii.read(path + 'SED_ECPL_pts.ecsv')
data_hap_ECPL = ascii.read(path + 'Spec_HYPEXP.ecsv')
# data_pa_PL    = ascii.read(path + 'SED_PL_pts.ecsv')
data_old      = ascii.read(path + 'spectrum_2005.ecsv')
data_pa_ECPL  = ascii.read(path + 'SED_STD_ECPL_pts.ecsv')
data_pa_PL    = ascii.read(path + 'SED_STD_PL_pts.ecsv')

# figure1 = naima.plot_data(data_pa_ECPL, e_unit=u.TeV, sed=False)
#
# naima.plot_data(data_hap_ECPL, e_unit=u.TeV, sed=False)
#
# naima.plot_data(data_old, e_unit=u.TeV, sed=False)
#
#
# figure3 = naima.plot_data(data_pa_ECPL, e_unit=u.TeV, sed=False)
# naima.plot_data(data_hap_ECPL, e_unit=u.TeV, sed=False, figure=figure3)
# naima.plot_data(data_old, e_unit=u.TeV, sed=False, figure=figure3)
#
#
# plt.show()


# model for particle distribution
ECPL     = naima.models.ExponentialCutoffPowerLaw( 1e36 * u.Unit('1/eV'), 1 * u.TeV, 2.1, 13 * u.TeV)

# Inverse Compton

IC       = naima.models.InverseCompton( ECPL, seed_photon_fields=['CMB', 'FIR', 'NIR'])
IC.particle_distribution.index = 2.2


# Synchrotron

SYN = naima.models.Synchrotron( ECPL, B =17 * u.uG)

# Emission spectra from the radiative models --> Flux and SED

spectrum_energy = np.logspace( -1, 14, 1000) * u.TeV
sed_IC  = IC.sed( spectrum_energy, distance = 4.0 * u.kpc)
sed_SYN = SYN.sed(spectrum_energy, distance = 4.0 * u.kpc)



#Model building

def lnprior(pars):
    lnprior = naima.uniform_prior(pars[0], 0., np.inf) \
            + naima.uniform_prior(pars[1], -1, 5)
    return lnprior

def model(pars, data):
    amplitude = pars[0] / u.TeV
    alpha = pars[1]
    e_cutoff = (10**pars[2]) * u.TeV

    ECPL = naima.models.ExponentialCutoffPowerLaw(amplitude, 10*u.TeV, alpha, e_cutoff)
    IC   = naima.models.InverseCompton(ECPL, seed_photon_fields=['CMB',
                        ['FIR', 45 * u.K, 0.42 * u.TeV / u.cm**3], 'NIR'])

    return IC.flux(data, distance=1.0*u.kpc)


# Initial Conditions
p0     = np.array((1e33, 2.0, np.log10(30)))
labels = ['norm', 'index', 'log10(cutoff)']

imf = naima.InteractiveModelFitter(model, p0, data=data_pa_ECPL, labels=labels, sed=False)

# posterior DF
sampler, pos = naima.run_sampler(data_table = data_pa_ECPL, p0=p0, model=model, prior=lnprior, nwalkers=128, nburn=50, nrun=10, threads=4)

naima.plot_fit(sampler, sed=False, e_unit=u.TeV )


naima.save_diagnostic_plots(path + 'MSH_15-52_PA_ECPL', sampler, blob_labels=['Spectrum', 'Electron energy distribution', '$W_e (E_e>1$ TeV)'])

naima.save_results_table(path + 'MSH15-52_naima_fit', sampler)





#Plot

plt.figure(figsize=(15,10))

plt.rc('font', family='sans')
plt.rc('mathtext', fontset='custom')

for seed, ls in zip(['CMB', 'FIR', 'NIR'], ['-','--',':']):
    sed = IC.sed(spectrum_energy, seed=seed, distance=1.5*u.kpc)
    plt.loglog(spectrum_energy,sed,lw=1,
            ls=ls,label='IC ({0})'.format(seed),c='0.25')


plt.loglog(spectrum_energy, sed_IC, lw=2, label='IC (total)', c=naima.plot.color_cycle[0])
#plt.loglog(spectrum_energy, sed_SYN, lw=2, label='Sync', c=naima.plot.color_cycle[1])
plt.xlabel('Photon energy [{0}]'.format(spectrum_energy.unit.to_string('latex_inline')))
#plt.ylabel('$E^2 dN/dE$ [{0}]'.format(sed_SYN.unit.to_string('latex_inline')))
plt.ylim(1e-13, 1e-6)
plt.tight_layout()
plt.legend(loc='lower left')

plt.show()


