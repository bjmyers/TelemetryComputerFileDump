# Updated to take into account April 2019 optic prescription changes.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib
import matplotlib.patches as patches
import sys
import pdb
from copy import deepcopy
from tqdm import tqdm
import astropy.units as u
import astropy.constants as c
from prettytable import PrettyTable
import soxs
from scipy import interpolate
import scipy.io as sio
import textwrap

sys.path.append('C:/python')
import PyXFocus.sources as sources
import PyXFocus.transformations as trans
import PyXFocus.surfaces as surfaces
import PyXFocus.analyses as analyses
import PyXFocus.conicsolve as conic
import OGRE.ogre_routines as ogre

# Configure plotting.
plt.style.use('seaborn')
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'


# Define default optic assembly misalignments.
oa_pitch = np.radians(120./3600)  # [rad]
oa_roll = np.radians(900./3600)  # [rad]
oa_yaw = np.radians(120./3600)  # [rad]
oa_x = 1.0  # [mm]
oa_y = 1.0  # [mm]
oa_z = 0.5  # [mm]
oa_misalign = [oa_pitch, oa_roll, oa_yaw, oa_x, oa_y, oa_z]

# Define center of optic assembly in z.
oa_z_cen = 3350.

# Define default jitter.
j = 0.42 * u.arcsec

# Create grating-to-grating misalignments.
g_pitch = np.radians(30./3600) / 1. * u.rad  # [rad]
g_roll = np.radians(15./3600) / 1. * u.rad  # [rad]
g_yaw = np.radians(30./3600) / 1. * u.rad  # [rad]
g_x = 0.127 * u.mm  # [mm]
g_y = 0.127 * u.mm  # [mm]
g_z = 0.127 * u.mm  # [mm]
g_misalign = ogre.create_g_misalign(g_pitch, g_roll, g_yaw, g_x, g_y, g_z)

# Create stack-to-stack misalignments.
s_pitch = np.radians(30./3600) * u.rad  # [rad]
s_roll = np.radians(20./3600) * u.rad  # [rad]
s_yaw = np.radians(60./3600) * u.rad  # [rad]
s_x = 0.254 * u.mm  # [mm]
s_y = 0.254 * u.mm  # [mm]
s_z = 0.127 * u.mm # [mm]
s_misalign = ogre.create_s_misalign(s_pitch, s_roll, s_yaw, s_x, s_y, s_z)

# Create module-to-array misalignments [pitch, roll, yaw, x, y, z].
m_misalign = np.array(
    [np.radians(120./3600), np.radians(120./3600), np.radians(120./3600),
     0.5, 0.5, 0.5])

det_z = 0.127


def interpolate_data(x_orig, y_orig, x_new, k=3, ext=3):
    '''
    Interpolates some data (x,y) over some new x-data.

    Parameters
    ----------
    x_orig : array-like
        Independent variable of data to be interpolated.
    y_orig : array-like
        Dependent variable of data to be interpolated.
    x_new : array-like
        New positions to interpolate data at.
    k : int, optional
        Degree of the smoothing spline. Must be 1 <= k <= 5.
    ext : int or str, optional
        Controls the extrapolation mode for elements not in the interval
        defined by the knot sequence.
            - If ext=0 or 'extrapolate', return the extrapolated value.
            - If ext=1 or 'zeros', return 0
            - If ext=2 or 'raise', raise a ValueError
            - If ext=3 of 'const', return the boundary value.
        The default value is 0.


    Returns
    -------
    y_new : array-like
        Interpolated values at specified positions.
    '''

    # Check to ensure data is the same length.
    if len(x_orig) != len(y_orig):
        raise ValueError('Provided data must have same length.')

    # Interpolate data.
    interp = interpolate.InterpolatedUnivariateSpline(
        x_orig, y_orig, k=k, ext=ext)

    return interp(x_new)


def generate_capella_spect(emin=0.2*u.keV, emax=1.25*u.keV, nbins=100000,
                           model='Audard2001'):
    '''
    Generates Capella spectrum from three-temperature thermal plasma
    model set forth by M. Audard et al. 2001.

    Parameters
    ----------
    emin : Quantity
        Minimum energy to model Capella spectrum.
    emax : Quantity
        Maximum energy to model Capella spectrum.
    nbins : int
        Number of energy points to sample between minimum (emin) and
        maximum energy (emax).

    Returns
    -------
    wv_bins : numpy.array
        Wavelength bins from the spectral simulation. Note: Bins are
        evenly distributed in energy space, but not in wavelength space.
        [Angstrom]
    flux : numpy.array
        Array giving the flux of the simulated Capella spectrum between
        two wavelength bins. [phot / s / cm2 / Angstrom]
    '''

    # Initiate ApecGenerator object.
    agen = soxs.ApecGenerator(
        emin.to('keV', equivalencies=u.spectral()),
        emax.to('keV', equivalencies=u.spectral()), nbins)

    if model is 'Audard2001':
        # Define plasma temperatures.
        t1 = 0.159  # [keV]
        t2 = 0.593  # [keV]
        t3 = 1.14  # [keV]

        # Define emission measures.
        em1 = 10**51.90  # [cm^-3]
        em2 = 10**52.86  # [cm^-3]
        em3 = 10**51.81  # [cm^-3]

        # Define other parameters.
        abund = 1.
        redshift = 0.
        dist = 13.159 * u.parsec
        dist = dist.to('cm').value
        norm1 = 1e-14 * em1 / 4 / np.pi / dist**2
        norm2 = 1e-14 * em2 / 4 / np.pi / dist**2
        norm3 = 1e-14 * em3 / 4 / np.pi / dist**2

        # Create spectra objects.
        spec_t1 = agen.get_spectrum(t1, abund, redshift, norm1)
        spec_t2 = agen.get_spectrum(t2, abund, redshift, norm2)
        spec_t3 = agen.get_spectrum(t3, abund, redshift, norm3)

        # Sum spectra.
        spec = spec_t1 + spec_t2 + spec_t3

        # Apply foreground extinction.
        nH = 0.00017988709151287835
        spec.apply_foreground_absorption(nH)

    elif model is 'Brickhouse2000':
        # Define plasma temperatures.
        t1 = (10**6.8 * u.Kelvin * c.k_B).to('keV')  # [keV]
        t2 = (10**7.0 * u.Kelvin * c.k_B).to('keV')  # [keV]

        # Define emission measures.
        em1 = 10**53  # [cm^-3]
        em2 = 10**52.7  # [cm^-3]

        # Define other parameters.
        abund = 1.
        redshift = 0.
        dist = 13.159 * u.parsec
        dist = dist.to('cm').value
        norm1 = 1e-14 * em1 / 4 / np.pi / dist**2
        norm2 = 1e-14 * em2 / 4 / np.pi / dist**2

        # Create spectra objects.
        spec_t1 = agen.get_spectrum(t1, abund, redshift, norm1)
        spec_t2 = agen.get_spectrum(t2, abund, redshift, norm2)

        # Sum spectra.
        spec = spec_t1 + spec_t2

        # Apply foreground extinction.
        nH = 0.00018
        spec.apply_foreground_absorption(nH)

    else:
        raise ValueError('Please provide a valid model name!')

    # Copy energy and flux arrays.
    e_bins = deepcopy(spec.ebins)
    e_flux = deepcopy(spec.flux)

    # Convert to wavelength.
    wave = deepcopy((e_bins).to('Angstrom', equivalencies=u.spectral()))

    # Convert flux from keV^-1 to Angstrom^-1
    flux = e_flux * np.mean(np.diff(e_bins))
    flux = flux / np.diff(wave) * -1

    wv_flux = deepcopy(flux)
    wv_bins = deepcopy(wave)

    return wv_bins, wv_flux


def get_optic_geom_area(spacer_width=4.*u.mm, azimuthal_gap=4*u.mm):
    '''
    Returns geometric area of OGRE optic.

    Parameters
    ----------
    spacer_width : Quantity or NoneType
        Defines bonding spacer width used for construct OGRE mirror
        module.
    azimuthal_gap : Quantity or NoneType
        Defines azimuthal gap between adjacent mirror segments in OGRE
        mirror module.

    Returns
    -------
    geom_area : float
        Geometric collecting area of the OGRE mirror module.
    '''

    # Catch if type of spacer_width is not what is expected.
    if type(spacer_width) is u.quantity.Quantity:
        pass
    elif spacer_width is None:
        pass
    else:
        raise TypeError(textwrap.fill(textwrap.dedent('''
            Spacer width must have an astropy.unit.quantity.Quantity object
            associated with it or it must be a NoneType.
            Please ensure that one of these two conditions have been
            met.''')))

    if type(azimuthal_gap) is u.quantity.Quantity:
        pass
    elif azimuthal_gap is None:
        pass
    else:
        raise TypeError(textwrap.fill(textwrap.dedent('''
            Azimuthal gap must have an astropy.unit.quantity.Quantity object
            associated with it or it must be a NoneType.
            Please ensure that one of these two conditions have been
            met.''')))

    # Define needed parameters and give them units.
    rp_front = ogre.rp_front
    rp_back = ogre.rp_back
    r_int = ogre.r_int

    # Define base geometric area or each shell.
    indiv_geom_area = (np.pi*rp_front**2 - np.pi*rp_back**2)

    if spacer_width is not None:
        # What fraction of the 360 degree is a single spacer taking up?
        ang_frac = np.arctan(
            (spacer_width / r_int).to('')) / (2 * np.pi * u.rad)
        # There will be 24 spacers per shell.
        tot_ang_frac = ang_frac * 24.
        # Calculate geometric collecting area of primary mirrors.
        indiv_geom_area *= (1 - tot_ang_frac)

    if azimuthal_gap is not None:
        # What fraction of the 360 degree is a single spacer taking up?
        ang_frac = np.arctan(
            (azimuthal_gap / r_int).to('')) / (2 * np.pi * u.rad)
        # There will be twelve gaps between mirrors in each shell.
        tot_ang_frac = ang_frac * 12.
        # Update geometric collecting area of primary mirrors.
        indiv_geom_area *= (1 - tot_ang_frac)

    # Sum geometric area.
    geom_area = np.sum(indiv_geom_area)

    return(geom_area.to('cm2'))


def get_primary_reflectivity(wave=None):
    '''
    Return primary reflectivity at specified wavelengths.

    Parameters
    ----------
    wave : None or Quantity
        Wavelengths to return reflectivity of the OGRE primary mirrors.
        If None, will return reflectivity values at the tabulated
        wavelengths.

    Returns
    -------
    wave : Quantity
        Wavelengths corresponding to the reflectivity values. If
        wavelength values are provided initially, this will be the same
        as those provided.
    reflectivity : Quantity
        Reflectivity values at each returned wavelength value.
    '''
    # Define filenames.
    prim_fname = 'C:/python/OGRE/Reflectivity/alpha-p-rev5-2nm-sio2-layer.dat'

    # Load in reflectivity data.
    prim_data = np.genfromtxt(
        prim_fname, skip_header=2, delimiter='  ').transpose()

    wv = deepcopy(prim_data[0]) * u.nm
    prim_reflect = u.Quantity(deepcopy(prim_data[1]))

    if wave is None:
        return wv, prim_reflect
    else:
        return wave, interpolate_data(
            wv.value, prim_reflect.value, wave.to(wv.unit)) * prim_reflect.unit


def get_secondary_reflectivity(wave=None):
    '''
    Return secondary reflectivity at specified wavelengths.

    Parameters
    ----------
    wave : None or Quantity
        Wavelengths to return reflectivity of the OGRE secondary mirrors.
        If None, will return reflectivity values at the tabulated
        wavelengths.

    Returns
    -------
    wave : Quantity
        Wavelengths corresponding to the reflectivity values. If
        wavelength values are provided initially, this will be the same
        as those provided.
    reflectivity : Quantity
        Reflectivity values at each returned wavelength value.
    '''
    # Define filenames.
    sec_fname = 'C:/python/OGRE/Reflectivity/alpha-h-rev5-2nm-sio2-layer.dat'

    # Load in reflectivity data.
    sec_data = np.genfromtxt(
        sec_fname, skip_header=2, delimiter='  ').transpose()

    wv = deepcopy(sec_data[0]) * u.nm
    sec_reflect = u.Quantity(deepcopy(sec_data[1]))

    if wave is None:
        return wv, sec_reflect
    else:
        return wave, interpolate_data(
            wv.value, sec_reflect.value, wave.to(wv.unit)) * sec_reflect.unit


def get_optic_eff_area(wave=None, spacer_width=4.*u.mm, azimuthal_gap=4*u.mm):
    '''
    Returns effective area for OGRE optic. Effective area is geometric
    area multiplied by the reflectivity of the primary mirrors and secondary
    mirrors.

    Here, we assume that we have a thick Si mirror (unpolarized
    incident photons), with a 2 nm thick Si02 layer
    (https://www.svmi.com/custom-film-coatings/thermal-oxide/). This
    function loads reflectivity data for primary
    and secondary mirrors of OGRE mirror module.

    Parameters
    ----------
    wave : None or Quantity
        Wavelength array over which to calculate effective area.
    spacer_width : Quantity or NoneType
        Defines bonding spacer width used for construct OGRE mirror
        module.
    azimuthal_gap : Quantity or NoneType
        Defines azimuthal gap between adjacent mirror segments in
        OGRE mirror module.

    Returns
    -------
    wave : Quantity
        Wavelength values corresponding to the calculated effective area.
    eff_area : Quantity
        Effective collecting area of the OGRE mirror module including
        reflectivity losses.
    '''


    # Get primary & secondary reflectivity.
    wv_prim, prim_reflect = get_primary_reflectivity(wave=wave)
    sec_reflect = get_secondary_reflectivity(wave=wave)[1]

    # Get geometric collecting area of OGRE optic.
    geom_area = get_optic_geom_area(spacer_width=spacer_width,
        azimuthal_gap=azimuthal_gap)

    if wave is None:
        # Calculate effective area.
        eff_area = geom_area * prim_reflect * sec_reflect
        return wv_prim, eff_area
    else:
        # Calculate effective area.
        eff_area = geom_area * prim_reflect * sec_reflect

        return wave, eff_area


def get_grating_throughput_frac(mid_width=5*u.mm, edge_support_width=5*u.mm,
                                grat_thickness=ogre.grating_thickness,
                                rib_width=2*u.mm, rib_num=3,
                                module_width=150*u.mm):
    '''
    Returns percentage of incident photons from optic that will interact
    with gratings in the OGRE grating array.

    Parameters
    ----------
    mid_width : Quantity
        Width between two grating stacks in grating module.
    edge_support_width : Quantity
        Width of left/right outside support in grating module.
    grat_thickness : Quantity
        Thickness of grating in OGRE grating stack.
    rib_width : Quantity
        Width of individual rib on wedged grating substrate.
    rib_num : Quantity
        Number of ribs on wedged grating substrate.
    module_width : Quantity
        Total width of OGRE grating module, centered at x=0.

    Returns
    -------
    throughput_frac : float
        Fraction of photons that interact with gratings in the OGRE
        grating array.
    '''

    # Create rays.
    rays = ogre.create_rays(n=1000000)

    # Pass them through the optic with Beckmann + Gaussian scatter.
    rays = ogre.ogre_mirror_module(rays, scatter='bg')

    # Move them to the mean grating center point.
    z_mean = (np.mean(ogre.zg_cen_left) + np.mean(ogre.zg_cen_right)) / 2
    trans.transform(rays, 0, 0, z_mean.to('mm').value, 0, 0, 0)
    surfaces.flat(rays)

    # Define starting number.
    start_num = len(rays[0])

    # Find photons that hit in between grating stacks.
    center_ind = np.where(np.abs(rays[1]) < mid_width.to('mm').value/2)[0]

    # Find photons with x-positions g.t. maximal grating stack x-extent.
    outside_ind = np.where(np.abs(rays[1]) > (
        module_width/2 - edge_support_width).to('mm').value)[0]

    # Find photons that hit ribs.
    grating_width = module_width/2 - edge_support_width - mid_width/2

    # Find rib positions.
    rib_cen_positions = np.linspace(
        rib_width.to('mm').value/2,
        (grating_width - rib_width/2).to('mm').value, rib_num) * u.mm + mid_width/2
    rib_ind = np.array([], dtype=int)
    for rc in rib_cen_positions:
        ri = np.where((np.abs(rays[1]) > (rc - rib_width/2).to('mm').value)
                      & (np.abs(rays[1]) < (rc + rib_width/2).to('mm').value))[0]
        rib_ind = np.concatenate((rib_ind, ri))

    # Account grating thickness.
    thick_ind = np.array([], dtype=int)
    for r in ogre.rg_cen.to('mm').value:
        # Find which rays will hit grating thickness.
        ind = np.where(
            (rays[2] < r) & (rays[2] > r - grat_thickness.to('mm').value))[0]
        # Add to array list.
        thick_ind = np.concatenate((thick_ind, ind))

    # Combine all rays that will be occulted.
    occ_ind = np.concatenate((center_ind, outside_ind, rib_ind, thick_ind))

    # Remove occulted rays.
    rays = [np.delete(r, occ_ind) for r in rays]

    # What is ending number?
    end_num = len(rays[1])

    # print('Starting #: ' + str(start_num))
    # print('Ending #: ' + str(end_num))

    # Fraction of photons making it through?
    frac = float(end_num) / start_num

    return frac


def get_corner_throughput_frac(module_width=150*u.mm):
    '''
    Returns percentage of incident photons from optic that will not
    interact with gratings in the OGRE grating array and will instead
    make it to the central CCD.

    Parameters
    ----------
    module_width : Quantity
        Total width of OGRE grating module, centered at x=0.

    Returns
    -------
    throughput_frac : float
        Fraction of photons that do not interact with the OGRE
        grating array.
    '''

    # Define rays.
    rays = ogre.create_rays()

    # Pass through optic.
    rays = ogre.ogre_mirror_module(rays, scatter='bg')

    # Move to center of OGRE grating module.
    z_cen = np.mean(np.concatenate((ogre.zg_cen_left, ogre.zg_cen_right)))
    trans.transform(rays, 0, 0, z_cen, 0, 0, 0)
    surfaces.flat(rays)

    # Find starting fraction.
    start_frac = len(rays[1])

    # How many get through.
    ind = np.where(abs(rays[1]) > module_width.to('mm').value / 2)[0]

    return len(ind) / float(start_frac)


def get_total_grating_efficiency(wave=None):
    '''
    Loads in up-to-date grating efficiency data and returns total grating
    diffraction efficiency.

    Currently, this data is from a PCGrate calculation of a 29.5 degree
    blazed grating at a graze angle of 1.5 degrees in the Littrow
    mounting. The grating is Au-coated and has a groove period of 160 nm.

    Parameters
    ----------
    wave : None or Quantity
        Array of wavelengths to return total grating efficiency over.
        If None, will return wavelength over range calculated by PCGrate.

    Returns
    -------
    wave : Quantity
        Array of wavelengths over which the efficiency data is
        calculated.
    eff : Quantity
        Array of total grating efficiency values.
    '''


    # Read in efficiency data.
    eff_data = np.genfromtxt(
        'C:/Users/Ben/Documents/OGRE/Efficiency/data_ogre_blazed_30deg.txt').transpose()

    # Find indices for given order.
    ind1 = np.where(eff_data[0] == -1.)
    ind2 = np.where(eff_data[0] == -2.)
    ind3 = np.where(eff_data[0] == -3.)
    ind4 = np.where(eff_data[0] == -4.)

    # Define wavelength for each order.
    wave1 = eff_data[1][ind1]
    wave2 = eff_data[1][ind2]
    wave3 = eff_data[1][ind3]
    wave4 = eff_data[1][ind4]

    # Find minimum and maximum wavelength values.
    wave_vals = np.concatenate([wave1, wave2, wave3, wave4])

    # Define wavelength array.
    wave_arr = np.unique(wave_vals) * u.nm

    # Calculate total diffraction efficiency in n=1,2,3,4.
    total_diff_eff = np.zeros(len(wave_arr), dtype=np.float64)
    diff_effs = np.concatenate([eff_data[2][ind1], eff_data[2][ind2],
                                eff_data[2][ind3], eff_data[2][ind4]])

    for i in range(len(wave_arr)):
        ind = np.where(wave_vals == wave_arr[i].value)[0]
        tot_eff = np.sum(diff_effs[ind])
        total_diff_eff[i] = tot_eff

    # Give tot_eff a unit.
    total_diff_eff = u.Quantity(total_diff_eff)

    if wave is not None:
        return wave, interpolate_data(
            wave_arr.value, total_diff_eff.value,
            wave.to(wave_arr.unit).value) * total_diff_eff.unit
    else:
        return wave_arr, total_diff_eff


def get_detector_filter_throughput(filt='poly-45nm-al-30nm-mesh', wave=None):
    '''
    Loads in throughput data for detector filter.

    Valid filters are "poly-45nm-al-30nm-mesh" & "al-50nm".

    Parameters
    ----------
    filter : str
        Name of filter data to load in.
    wave : None or Quantity
        Array of wavelengths to interpolate filter throughput over. If
        None, will return wavelength over range in original filter data.

    Returns
    -------
    wave : Quantity
        Array of wavelengths over which specified detector filter
        throughput is returned.
    filter_tput : Quantity
        Values of specified detector filter throughput at each wavelength
        value.
    '''

    if filt == 'poly-45nm-al-30nm-mesh':
        # Load in filter throughput.
        filter_data = sio.loadmat(
            'C:/Users/Ben/Documents/OGRE/Filters/APRA_filters.mat')
        filter_tput = filter_data['filter_5'].transpose()[0]

    elif filt =='al-50nm':
        # Load in filter throughput.
        filter_data = sio.loadmat(
            'C:/Users/Ben/Documents/OGRE/Filters/Al_filters.mat')
        filter_tput = filter_data['filter_1'].transpose()[0]

    else:
        # Raise value error if *filt* is not valid string.
        raise ValueError(textwrap.fill(textwrap.dedent('''
            Specified filter name was not expected. Please choose only a valid
            filter option from those listed in docstring.''')))

    # Load in energy data.
    energy_data = sio.loadmat('C:/Users/Ben/Documents/OGRE/Filters/energy.mat')
    energy = energy_data['energy'][0] * u.keV
    wave_arr = energy.to('nm', equivalencies=u.spectral())

    # Limit data b/c original data is a BEAST!
    ind_filt = np.where((wave_arr.value > 1.) & (wave_arr.value < 10.))
    filter_tput = deepcopy(filter_tput[ind_filt])
    wave_arr = deepcopy(wave_arr[ind_filt])

    # Give filter throughput a unit.
    filter_tput = u.Quantity(filter_tput)

    if wave is not None:
        return wave, interpolate_data(
            wave_arr[::-1].value, filter_tput[::-1].value,
            wave.to(wave_arr.unit).value, k=1) * filter_tput.unit

    else:
        return wave_arr, filter_tput


def get_detector_qe(wave=None):
    '''
    Returns quantum efficiency of OGRE EM-CCD.

    Parameters
    ----------
    wave : None or Quantity
        Wavelength values of which to interpolate detector quantum efficiency.

    Returns
    -------
    wave : Quantity
        Wavelength values corresponding to the returned detector QE values.
    det_qe : Quantity
        QE of OGRE EM-CCD at each wavelength value.
    '''

    # Load in QE data.
    qe_data = sio.loadmat(
        'C:/Users/Ben/Documents/OGRE/Filters/back_illuminated_CCD.mat')
    qe = qe_data['totalQEb2'][0] * u.count / u.photon
    energy_qe = qe_data['energy'][0] * u.keV

    # Filter only needed data.
    wv_qe = energy_qe.to('nm', equivalencies=u.spectral())
    ind_qe = np.where((wv_qe.value > 0.5) & (wv_qe.value < 10.))[0]
    wv_qe = deepcopy(wv_qe[ind_qe])
    qe = deepcopy(qe[ind_qe])

    if wave is not None:
        return wave, interpolate_data(
            wv_qe[::-1].value, qe[::-1].value,
            wave.to(wv_qe.unit).value, k=1) * qe.unit

    else:
        return wv_qe[::-1], qe[::-1]


def wave_values_from_capella_spect(rays, wave_start, wave_end):
    '''
    Creates an array of wavelengths and diffraction orders from
    SOXS-simulated Capella spectrum and grating diffraction efficiencies.
    The length of each array will be the same length as number of photons
    in the rays variable.

    Parameters
    ---------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        Rays object.
    wave_start : Quantity
        Lower wavelength in wavelength range to simulate.
    wave_end : Quantity
        Upper wavelength in wavelength range to simulate.

    Returns
    -------
    waves : Quantity
        List of wavelengths with the same length as len(rays[0]).
        Wavelengths will be representative of observed Capella spectrum.
    orders : numpy.ndarray
        Array of diffraction orders corresponding to the wavelengths also
        returned. Orders are chosen based on diffraction efficiencies.
    '''
    # Find how many rays we'll be analyzing.
    n_rays = len(rays[0])

    # Generate Capella spectrum.
    wv_bins, wv_flux = generate_capella_spect(emin=wave_end, emax=wave_start)

    # Redefine wv_flux as prob_vals.
    prob_vals = deepcopy(wv_flux)

    # Scale so that sum is unity.
    prob_vals /= np.sum(prob_vals)

    # Move from wavelength bin edges to wavelength bin centers.
    bin_diff = np.diff(wv_bins)
    wv_cen = wv_bins[:-1] + bin_diff/2  # [Angstrom]

    # Randomly draw wavelengths from prob_vals distribution.
    waves = np.random.choice(
        wv_cen, n_rays, replace=True, p=prob_vals) * u.Angstrom

    # Read in efficiency data.
    eff_data_fname = 'C:/Users/Ben/Documents/OGRE/Efficiency/data_ogre_blazed_30deg.txt'
    eff_data = np.genfromtxt(eff_data_fname).transpose()

    # Find indices which correspond to first through fourth orders.
    first_ind = np.where(eff_data[0] == -1.)[0]
    second_ind = np.where(eff_data[0] == -2.)[0]
    third_ind = np.where(eff_data[0] == -3.)[0]
    fourth_ind = np.where(eff_data[0] == -4.)[0]

    # Fit each of these orders to a spline.
    # Wavelength values are [nm].
    # ext=1 returns 0 if outside of provided range.
    first_interp = interpolate.InterpolatedUnivariateSpline(
        eff_data[1][first_ind], eff_data[2][first_ind], k=2, ext=1)
    second_interp = interpolate.InterpolatedUnivariateSpline(
        eff_data[1][second_ind], eff_data[2][second_ind], k=2, ext=1)
    third_interp = interpolate.InterpolatedUnivariateSpline(
        eff_data[1][third_ind], eff_data[2][third_ind], k=2, ext=1)
    fourth_interp = interpolate.InterpolatedUnivariateSpline(
        eff_data[1][fourth_ind], eff_data[2][fourth_ind], k=2, ext=1)

    # Assign order to each wavelength based on probabilities.
    orders = np.arange(1, 5)
    order_vals = []
    for wave in waves.to('nm').value:  # Convert to [nm] for interp.
        prob = np.array([first_interp(wave), second_interp(wave),
                         third_interp(wave), fourth_interp(wave)])
        prob = prob / np.sum(prob)  # Make probs sum to 1.
        order_vals.append(
            np.random.choice(orders, 1, replace=True, p=prob)[0])
    order_vals = np.array(order_vals)

    orders = deepcopy(order_vals)

    return waves, orders


def get_central_ccd_effective_area(wave, det_filt='poly-45nm-al-30nm-mesh',
                                   module_width=150*u.mm, spacer_width=4*u.mm,
                                   azimuthal_gap=4*u.mm):
    '''
    Calculates the effective area of optic on central CCD for given
    wavelength values.

    Parameters
    ----------
    wave : Quantity
        Wavelength values over which to calculate effective area.
    det_filt : str
        Filter throughput values to use. Default is
        'poly-45nm-al-30nm-mesh'. Could also specify 'al-50nm'.
    spacer_width : Quantity or NoneType
        Defines bonding spacer width used to construct OGRE mirror
        module.
    azimuthal_gap : Quantity or NoneType
        Defines azimuthal gap between adjacent mirror segments in
        OGRE mirror module.

    Returns
    -------
    wave : Quantity
        Wavelengths corresponding to the returned effective area values.
    eff_area : Quanity
        Effective area at corresponding wavelength values for the OGRE
        spectrometer.
    '''

    if wave is None:
        raise ValueError('Must specify a valid wavelength array!')

    # Calculate total effective area.
    eff_area = get_optic_eff_area(
        wave=wave, spacer_width=spacer_width, azimuthal_gap=azimuthal_gap)[1]
    eff_area *= get_corner_throughput_frac(module_width=module_width)
    eff_area *= get_detector_filter_throughput(wave=wave, filt=det_filt)[1]

    return wave, eff_area


def get_spectrum_on_central_ccd(wave_start=2.0*u.keV, wave_end=0.2*u.keV,
                                obs_time=300*u.second):
    '''
    Returns the spectrum of photons that do not interact with the OGRE
    grating module and pass to central detector.

    Parameters
    ----------
    wave_start : Quantity
        Lower wavelength in wavelength range to simulate.
    wave_end : Quantity
        Upper wavelength in wavelength range to simulate.
    '''

    # Get Capella spectrum.
    wv_bins, wv_flux = generate_capella_spect(emin=wave_end, emax=wave_start)

    # Move from wavelength bin edges to wavelength bin centers.
    bin_diff = np.diff(wv_bins)
    wv_cen = wv_bins[:-1] + bin_diff/2  # [Angstrom]

    # Get effective area.
    eff_area = get_central_ccd_effective_area(wv_cen)[1]

    # Multiply by effective area.
    wv_flux *= eff_area

    # Multiply by observation time.
    wv_flux *= obs_time

    # Get detector qe.
    det_qe = get_detector_qe(wave=wv_cen)[1]
    wv_flux *= det_qe

    return wv_bins, wv_flux




def get_effective_area(wave, det_filt='poly-45nm-al-30nm-mesh',
                       mid_width=5*u.mm, edge_support_width=5*u.mm,
                       grat_thickness=ogre.grating_thickness,
                       rib_width=2*u.mm, rib_num=3, module_width=150*u.mm,
                       spacer_width=4.*u.mm, azimuthal_gap=4*u.mm):
    '''
    Calculates the effective area for given wavelength values.

    Parameters
    ----------
    wave : Quantity
        Wavelength values over which to calculate effective area.
    det_filt : str
        Filter throughput values to use. Default is
        'poly-45nm-al-30nm-mesh'. Could also specify 'al-50nm'.
    mid_width : Quantity
        Width between two grating stacks in grating module.
    edge_support_width : Quantity
        Width of left/right outside support in grating module.
    grat_thickness : Quantity
        Thickness of grating in OGRE grating stack.
    rib_width : Quantity
        Width of individual rib on wedged grating substrate.
    rib_num : int
        Number of ribs on wedged grating substrate.
    module_width : Quantity
        Total width of OGRE grating module, centered at x=0.
    spacer_width : Quantity or NoneType
        Defines bonding spacer width used to construct OGRE mirror
        module.
    azimuthal_gap : Quantity or NoneType
        Defines azimuthal gap between adjacent mirror segments in
        OGRE mirror module.

    Returns
    -------
    wave : Quantity
        Wavelengths corresponding to the returned effective area values.
    eff_area : Quanity
        Effective area at corresponding wavelength values for the OGRE
        spectrometer.
    '''

    if wave is None:
        raise ValueError('Must specify a valid wavelength array!')

    # Calculate total effective area.
    eff_area = get_optic_eff_area(
        wave=wave, spacer_width=spacer_width, azimuthal_gap=azimuthal_gap)[1]
    eff_area *= get_grating_throughput_frac(
        mid_width=mid_width, edge_support_width=edge_support_width,
        grat_thickness=grat_thickness, rib_width=rib_width, rib_num=rib_num,
        module_width=module_width)
    eff_area *= get_total_grating_efficiency(wave=wave)[1]
    eff_area *= get_detector_filter_throughput(wave=wave, filt=det_filt)[1]

    return wave, eff_area


def get_total_counts(obs_time=270*u.second, wave_start=10*u.Angstrom,
                     wave_end=20*u.Angstrom):
    '''
    Gets total estimated counts observed by OGRE spectrometer in
    specified bandpass.

    Parameters
    ----------
    obs_time : Quantity
        Observation time.
    wave_start : Quantity
        One-half of bandpass of interest.
    wave_end : Quantity
        Other half of bandpass of interest.

    Returns
    -------
    tot_counts : Quantity
        Estimated count rate between wave_start and wave_end.
    '''

    # Get Capella spectrum.
    wv_bins, wv_flux = generate_capella_spect(
        emin=wave_end, emax=wave_start)

    # Move from wavelength bin edges to wavelength bin centers.
    bin_diff = np.diff(wv_bins)
    wv_cen = wv_bins[:-1] + bin_diff/2  # [Angstrom]

    # Caculate effective area of OGRE spectrometer.
    eff_area = get_effective_area(wv_cen)[1]

    # Get observed flux.
    obs_flux = wv_flux * eff_area * get_detector_qe(wave=wv_cen)[1]

    # Account for observation time.
    obs_spec = obs_flux * obs_time

    # Multiply by bin width.
    tot_counts = obs_spec * -np.diff(wv_bins)

    return np.sum(tot_counts)


def ogre_simulation_one_channel(waves, orders, j=j, oa_misalign=oa_misalign,
                               g_misalign=g_misalign, s_misalign=s_misalign,
                               m_misalign=m_misalign, cen_width=5*u.mm,
                               rib_width=2*u.mm, rib_num=3,
                               edge_width=5.*u.mm, mod_width=150*u.mm):
    '''
    Simulate rays going through the OGRE spectrometer. This function will
    just simulate rays passing through the spectrometer with
    misalignments, but will not take throughput into consideration. A
    seperate function is needed to correct this output for throughput.

    Parameters
    ----------
    waves : Quantity
        Array of wavelengths to simulate. Each wavelength in the array
        will be associated with one ray (100 wavelengths = 100 rays).
    orders : numpy.ndarray
        Diffraction orders associated with wavelengths passed.
    j : Quantity
        FWHM of the jitter throughout the observation.
    oa_misalign : 
    g_misalign : 
    s_misalign : 
    m_misalign : 

    Returns
    -------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        Rays on the focal plane of the OGRE spectrometer.
    waves : Quantity
        Wavelengths corresponding to each simulated ray.
    orders : array
        Diffracted order of each simulated ray.

    '''

    # Define initial rays.
    n = len(waves)
    rays = ogre.create_rays(10 * n)

    # Translate rays to center of optic assembly.
    rays[3] = np.ones(len(rays[1])) * oa_z_cen
    surfaces.flat(rays)

    # Misalign rays relative to optic assembly.
    trans.transform(rays, 0, 0, 0, oa_pitch, oa_yaw, 0)

    # Move rays to R_INT position. This is just a global offset.
    rays[3] = np.ones(len(rays[1])) * 3500.

    # Add jitter.
    xj = np.random.normal(0, j.to('rad').value, size=len(rays[4]))
    yj = np.random.normal(0, j.to('rad').value, size=len(rays[4]))

    rays[4] += xj
    rays[5] += yj
    rays[6] = np.sqrt(1 - rays[4]**2 - rays[5]**2)

    # Pass rays through optic.
    rays = ogre.ogre_mirror_module(rays, scatter='bg')

    # Randomly select indicies to keep based on count rate analysis.
    n_rays = len(rays[0])
    rand_inds = np.random.choice(
        np.arange(n_rays, dtype=int), len(waves), replace=False)

    # Keep only those selected indicies.
    rays = [r[rand_inds] for r in rays]

    # Add random period error effect.
    rand_fact = np.random.normal(1, 0.000094, size=len(waves))

    # Pass through grating moudle.
    rays, diff_ind = ogre.ogre_grating_module(
        rays, waves=waves*orders/rand_fact, g_misalign=g_misalign,
        s_misalign=s_misalign, m_misalign=m_misalign, return_diff_ind=True)

    # Keep only wavelength/orders that make it through grating module.
    waves = waves[diff_ind]
    orders = orders[diff_ind]

    # Occult rays that will hit ribs.
    rays, waves, orders = ogre.grating_rib_occultation(
        rays, waves=waves, orders=orders)

    # Transform coordinate system to be at center of rotation.
    trans.transform(rays, 0, 0, oa_z_cen, 0, 0, 0)

    # Propagate rays to misaligned plane, simulating center plane of optic assembly.
    surfaces.flat(rays)

    # Transform the rays back by the pitch angle to the optical axis coordinate system.
    trans.itransform(rays, 0, 0, 0, oa_pitch, oa_yaw, 0)

    # Put the rays back at the intersection plane of the ideal optic.
    rays[3] = np.ones(len(rays[3])) * oa_z_cen

    # Misalign in z-hat.
    trans.transform(rays, 0, 0, oa_z, 0, 0, 0)

    # Propagate to the focal plane.
    rays = ogre.ogre_focal_plane(
        rays, det_misalignments=[0, 0, 0, 0, 0, -det_z])

    return rays, waves, orders


def optic_reflectivity_correction(rays, waves, orders):
    '''
    Takes in a list of rays and removes some of the rays based upon
    reflectivity of primary and secondary mirrors.

    Parameters
    ----------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        PyXFocus rays variable.
    waves : Quantity
        Wavelengths corresponding to each ray in rays parameter.
    orders : numpy.ndarray
        Diffraction orders corresponding to each wavelength value.

    Returns
    -------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        PyXFocus rays variable after reflectivity correction.
    waves : Quantity
        Wavelengths corresponding to each ray in rays parameter after
        reflectivity correction.
    orders : numpy.ndarray
        Diffraction orders corresponding to each wavelength value after
        reflectivity correction.
    '''

    # Get total reflectivity off primary + secondary.
    r_prob = get_primary_reflectivity(wave=waves)[1].value * \
        get_secondary_reflectivity(wave=waves)[1].value

    # Step through array to determine which are reflected.
    r_inds = []
    for i in xrange(len(waves)):
        refl = -1
        # Randomly select if it is reflected. 0 = absorbed, 1 = refl.
        refl = int(np.random.choice([0,1], p=[1-r_prob[i], r_prob[i]]))
        if refl == 1:
            r_inds.append(i)

    # Keep only reflected rays.
    rays = [r[r_inds] for r in deepcopy(rays)]

    return rays, waves[r_inds], orders[r_inds]


def tot_diff_eff_correction(rays, waves, orders):
    '''
    Takes in a list of rays and removes some of the rays based upon
    the total diffraction efficiency of the gratings.

    Note: The relative populations in the different diffraction orders
    have already been set by a function such as
    wave_values_from_capella_spect().

    Parameters
    ----------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        PyXFocus rays variable.
    waves : Quantity
        Wavelengths corresponding to each ray in rays parameter.
    orders : numpy.ndarray
        Diffraction orders corresponding to each wavelength value.

    Returns
    -------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        PyXFocus rays variable after the total diffraction efficiency
        correction.
    waves : Quantity
        Wavelengths corresponding to each ray in rays parameter after
        the total diffraction efficiency correction.
    orders : numpy.ndarray
        Diffraction orders corresponding to each wavelength value after
        the total diffraction efficiency correction.
    '''

    # Get diffraction efficiency probability for wavelength values.
    d_prob = get_total_grating_efficiency(wave=waves)[1].value

    # Step through array to determine which are diffracted.
    d_inds = []
    for i in xrange(len(waves)):
        diff = -1
        # Randomly select if it is reflected. 0 = absorbed, 1 = refl.
        diff = int(np.random.choice([0, 1], p=[1-d_prob[i], d_prob[i]]))
        if diff == 1:
            d_inds.append(i)

    # Keep only reflected rays.
    rays = [r[d_inds] for r in deepcopy(rays)]

    return rays, waves[d_inds], orders[d_inds]


def det_filt_correction(rays, waves, orders):
    '''
    Takes in a list of rays and removes some of the rays based upon
    the detector filter transmission.

    Parameters
    ----------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        PyXFocus rays variable.
    waves : Quantity
        Wavelengths corresponding to each ray in rays parameter.
    orders : numpy.ndarray
        Diffraction orders corresponding to each wavelength value.

    Returns
    -------
    rays : [opd, x, y, z, l, m, n, nx, ny, nz]
        PyXFocus rays variable after the detector filter transmission
        correction.
    waves : Quantity
        Wavelengths corresponding to each ray in rays parameter after
        the detector filter transmission correction.
    orders : numpy.ndarray
        Diffraction orders corresponding to each wavelength value after
        the detector filter transmission correction.
    '''

    # Get diffraction efficiency probability for wavelength values.
    t_prob = get_detector_filter_throughput(wave=waves)[1].value

    # Step through array to determine which are diffracted.
    t_inds = []
    for i in xrange(len(waves)):
        trans = -1
        # Randomly select if it is reflected. 0 = absorbed, 1 = refl.
        trans = int(np.random.choice([0, 1], p=[1-t_prob[i], t_prob[i]]))
        if trans == 1:
            t_inds.append(i)

    # Keep only reflected rays.
    rays = [r[t_inds] for r in deepcopy(rays)]

    return rays, waves[t_inds], orders[t_inds]


def det_qe_correction(rays, waves, orders):
    '''
    Correct list of rays to take into account detector QE.
    '''
    # Get count probability for wavelength values.
    c_prob = get_detector_filter_throughput(wave=waves)[1].value

    # Step through array to determine which are diffracted.
    c_inds = []
    for i in xrange(len(waves)):
        count = -1
        # Randomly select if it is reflected. 0 = absorbed, 1 = refl.
        count = int(np.random.choice([0, 1], p=[1-c_prob[i], c_prob[i]]))
        if count == 1:
            c_inds.append(i)

    # Keep only reflected rays.
    rays = [r[c_inds] for r in deepcopy(rays)]

    return rays, waves[c_inds], orders[c_inds]


def simulate_ogre_observed_capella_spec(wave_start, wave_end, obs_time):
    '''
    Simulates an OGRE observed Capella spectrum.

    '''
    if type(wave_start) is not u.quantity.Quantity:
        raise ValueError()

    if type(wave_end) is not u.quantity.Quantity:
        raise ValueError()

    if type(obs_time) is not u.quantity.Quantity:
        raise ValueError()

    # Get Capella spectrum.
    wv_bins, wv_flux = generate_capella_spect(emin=wave_end, emax=wave_start)

    # Figure out how many photons are hitting optic aperture.
    optic_area = np.pi*(ogre.max_radial_extent**2 - ogre.min_radial_extent**2)
    int_flux = np.diff(wv_bins) * -1 * wv_flux * obs_time * optic_area
    int_flux = int_flux.to('ph').sum()

    # Find total flux.
    n_tot = int(int_flux.value)

    # Create wavelengths and diffraction orders.
    waves, orders = wave_values_from_capella_spect(
        [np.ones(n_tot)], wave_start, wave_end)

    # Pass through OGRE spectrometer.
    rays, waves, orders = ogre_simulation_one_channel(waves, orders)

    # Account for optic reflectivity.
    rays, waves, orders = optic_reflectivity_correction(rays, waves, orders)

    # Account for detector filter.
    rays, waves, orders = det_filt_correction(rays, waves, orders)

    # Account for total diffraction efficiency.
    rays, waves, orders = tot_diff_eff_correction(rays, waves, orders)

    # Account for detector filter throughput.
    rays, waves, orders = det_filt_correction(rays, waves, orders)

    # Account for detector QE.
    rays, waves, orders = det_qe_correction(rays, waves, orders)

    return rays, waves, orders
