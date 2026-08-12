import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.special import erf
from mendeleev import element
from tqdm import tqdm
from multiprocessing import Pool

from flux_history import FluxHistory

# --- Physical Constants ---
PROTON_MASS_MEV = 938.3
NEUTRON_MASS_MEV = 939.6
AVOGADRO_NUMBER = 6.022e23
MYR_PER_SECOND = 1/ (60 * 60 * 24 * 365 * 1e6)

# --- Binning setup ---
RECOIL_N_BINS= 201
RECOIL_ER_MIN_LOG_MEV= -2 
RECOIL_ER_MAX_LOG_MEV= 6
RECOIL_ENERGY_BINS_MEV = np.logspace(RECOIL_ER_MIN_LOG_MEV, RECOIL_ER_MAX_LOG_MEV, RECOIL_N_BINS)

LENGTH_N_BINS = 1000
LENGTH_MIN_LOG_NM = 1.5
LENGTH_MAX_LOG_NM = 5.5
TRACK_LENGTH_BINS_NM = np.logspace(LENGTH_MIN_LOG_NM, LENGTH_MAX_LOG_NM, LENGTH_N_BINS)

TYPICAL_DEPTH_MM = 0.001

# --- Utility Functions ---
def log_interp1d(xx, yy, kind='linear'):
    """
    Performs an interpolation in log-log space for better accuracy with wide-ranging data.

    Args:
        xx (np.ndarray): The x-coordinates of the data points.
        yy (np.ndarray): The y-coordinates of the data points.
        kind (str, optional): Specifies the kind of interpolation as a string. Defaults to 'linear'.

    Returns:
        function: An interpolation function that operates in log-log space.
    """
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interp1d(logx, logy, kind=kind, fill_value='extrapolate')
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

def calibrate_spectrum(x_bins, counts, x_scale_factor=1.0, y_scale_factor=1.0):
    """
    Applies an absolute calibration shift to the x-axis (track length) and y-axis (counts) of a spectrum.

    This function scales the track lengths by the given factor and then re-bins
    the counts onto the original binning structure using linear interpolation.

    Args:
        x_bins (np.ndarray): The bin edges for track length spectrum [nm].
        counts (np.ndarray): The array of track counts in each bin.
        x_scale_factor (float): The multiplicative factor to apply to the x-axis.
        y_scale_factor (float): The multiplicative factor to apply to the y-axis.

    Returns:
        np.ndarray: The calibrated track counts, re-binned onto the original
                    x_bins structure.
    """

    x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0

    x_scaled = x_mids * x_scale_factor

    calibrated_y = np.interp(x_mids, x_scaled, counts)
    calibrated_y *= (np.sum(counts) / np.sum(calibrated_y))

    return calibrated_y*y_scale_factor

def slice_spectrum(x_bins, counts, angular_pdf=None, phi_cut_deg=0., pit_width=500., bulk_etching_depth=TYPICAL_DEPTH_MM*1.e6, f_phi= lambda phi: 1., correction=True):
    """
    Applies Monte Carlo simulation of track slicing, accounting for geometrical, angular, 
    and experimental filtering effects (min/max measurable length).

    Args:
        x_bins (np.ndarray): The bin edges for the true track length spectrum R [nm].
        counts (np.ndarray): The array of true track counts N(R) in each bin.
        angular_pdf (np.ndarray, optional): Normalized array P(phi) for the angle distribution 
                                            of tracks relative to the surface normal. Defaults to isotropic (sin(phi)).
        phi_cut_deg (float): Angular filter threshold (tracks with phi < phi_cut are rejected). 
                                Set to 0.0 for highly-faithful plasma etching.
        pit_width (float): Typical width of the etched pit [nm]. 
                            Tracks with parallel footprint smaller than this threshold will be measured by this.
        bulk_etching_depth (float): Vertical development of the etching [nm]. 
        f_phi (callable): Function applied to the segment length L_seg * f_phi(phi). Corrects 
                            for anisotropic enlargement (e.g., set to lambda phi: 1.0 for plasma etching).
        correction (bool): If True, applies pit_width correction. If not, the pit_width correction is ignored.
    Returns:
        np.ndarray: The resulting measured track count histogram N(L_meas), normalized to 
                    the total input counts.
    """

    x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
    phi_cut_rad = np.deg2rad(phi_cut_deg)

    factor = int(x_bins[-1]/bulk_etching_depth)

    stat_factor = np.sum(counts) * (factor + 1)


    n_samples = stat_factor * 10**np.rint(8 - np.log10(stat_factor))

    samples = np.random.choice(x_mids, size=int(n_samples), p=counts/np.sum(counts))

    if angular_pdf:
        phi_grid = np.linspace(0, np.pi / 2, 1000)
        sampled_angles = np.random.choice(phi_grid, size=int(n_samples), p=angular_pdf / np.sum(angular_pdf))
    else:
        sampled_angles = np.random.uniform(low=0, high=np.pi / 2, size=int(n_samples))

    is_retained = sampled_angles >= phi_cut_rad

    samples_retained = samples[is_retained]
    phi_retained = sampled_angles[is_retained]

    sim_start_point = np.random.uniform(low = -factor*bulk_etching_depth, high = bulk_etching_depth, size=len(samples_retained))

    sim_end_point = sim_start_point + (samples_retained * np.sin(phi_retained))

    valid = (sim_end_point > bulk_etching_depth)

    angles_valid = phi_retained[valid]
    cut_sim_true_depth = sim_end_point[valid]

    depth = (cut_sim_true_depth - bulk_etching_depth)/np.sin(angles_valid)

    measured_samples = depth*f_phi(angles_valid)

    if correction:
        corrected_measurable_samples = np.where(measured_samples * np.cos(angles_valid) >= pit_width/2., (pit_width/2.)+measured_samples * np.cos(angles_valid), pit_width)
    else:
        corrected_measurable_samples = measured_samples

    hist_measured, _ = np.histogram(corrected_measurable_samples, bins=x_bins, density=False)

    hist_norm = hist_measured / 10**np.rint(8 - np.log10(stat_factor))

    return hist_norm

def detection_model_efficiency(x_bins, counts, precision, recall, model_mean, sigma_left, sigma_right=None, meas_error=1000.):
    """
    Multiplies counts (sliced) with the efficiency function for the detection model
    
    Args:
        x_bins (np.ndarray): The bin edges for the track length spectrum R [nm].
        counts (np.ndarray): The array of track counts N(R) in each bin.
        precision (float): Precision of the counting model.
        recall (float): Recall of the counting model.
        model_mean (float): Peak length for the efficiency distribution.
        sigma_left (float): Left standard deviation of the efficiency distirbution assuming asymmetric normal shape.
        sigma_right (float, optional): Right standard deviation of the efficiency distirbution assuming asymmetric normal shape.
        
    Returns:
        np.ndarray: Efficiency-corrected counts.
    """

    x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0

    eff = np.zeros_like(x_mids)

    if sigma_right is None:
        sigma_right = sigma_left

    center = np.argmin(np.abs(model_mean - x_mids))

    eff[:center] = np.exp(-(x_mids[:center] - model_mean)**2 / (2 * sigma_left**2))
    eff[center] = 1.0
    eff[center+1:] = np.exp(-(x_mids[center+1:] - model_mean)**2 / (2 * sigma_right**2))

    counts_with_efficiency = counts * eff * recall / precision

    counts_with_measure = smear_spectrum(counts_with_efficiency, len(x_bins)//2*2-1, meas_error/np.diff(x_bins)[0], meas_error/np.diff(x_bins)[0])

    return counts_with_measure

# --- Main Paleodetector Class ---
class Paleodetector:
    """
    A class to handle all physics calculations for a specific self.
    
    This class is initialized with a configuration dictionary and provides
    methods to calculate various signal and background components for
    paleo-detector analysis.
    """
    
    def __init__(self, mineral_config, total_age_kyr, overburden_history=None, flux_history=None, data_path_prefix="Data"):
        """
        Initializes the Paleodetector object and sets up its properties and data caches.

        Args:
            mineral_config (dict): A dictionary containing all properties of the mineral,
                                   such as name, composition, and nuclear data.
            total_age_kyr (float): The total age of the sample in thousands of years.
            overburden_history (dict, optional): A dictionary containing the overburden history data.
            flux_history (FluxHistory, optional): A FluxHistory instance describing the cosmic-ray
                                              primary/muon/neutron flux over time. Its own timeline
                                              runs 0 (present) -> negative (past); local t_kyr=0 here
                                              (start of data taking) corresponds to
                                              t_kyr=-total_age_kyr on that timeline, and local
                                              t_kyr=total_age_kyr (present) corresponds to t_kyr=0
                                              there. If None, a baseline-only FluxHistory is built
                                              (or loaded from cache) the first time it's needed.
            data_path_prefix (str, optional): The relative path to the main data directory.
                                              Defaults to "Data/".
        """
        self.config = mineral_config
        self.name = mineral_config['name']
        self.shortname = mineral_config["shortname"]
        self.composition = mineral_config['composition']
        self.data_path = data_path_prefix

        self.verbose = 1

        self.radiogenic_spectrum = self._radiogenic_spectrum()
        
        self.alpha_n_spectrum = self._alpha_n_spectrum()

        self.total_age_kyr = total_age_kyr

        self._overburden_interpolator = self._interpolate_overburden_history(overburden_history)

        self.flux_history = flux_history

        self._srim_cache = {}
        self._nuclear_data_cache = {}
        self._neutron_bkg_cache = {}
        self._depth_interpolators = {}
        self._secondary_n_spectrum = {}
        
        if self.verbose>0:
            print(f"Initialized Paleodetector: {self.name}")

    def set_flux_history(self, FluxHistory):
        self.flux_history = FluxHistory

    def _interpolate_overburden_history(self, overburden_history=None):
        """
        Interpolates the overburden history data to create a function that returns
        the overburden in meters water equivalent (mwe) for any given time.

        Args:
            overburden_history (dict, optional): A dictionary containing the overburden history data.
                If None, a default constant overburden of 0 mwe is used.
        """
        if isinstance(overburden_history, (int, float)):
            const = float(overburden_history)
            return interp1d([0.0, self.total_age_kyr], [const, const],
                             bounds_error=False, fill_value=const)
        elif overburden_history is None:
            return interp1d([0.0, self.total_age_kyr], [0.0, 0.0],
                             bounds_error=False, fill_value=0.0)

        initial_depth = overburden_history.get("initial_depth", 0.3)
        initial_density_g_cm3 = overburden_history.get("initial_density_g_cm3", 1.0)

        start_time_continuous_kyr = np.atleast_1d(overburden_history["start_time_continuous_kyr"]) if "start_time_continuous_kyr" in overburden_history else [0.0]
        end_time_continuous_kyr = np.atleast_1d(overburden_history["end_time_continuous_kyr"]) if "end_time_continuous_kyr" in overburden_history else [self.total_age_kyr]
        rate_continuous = np.atleast_1d(overburden_history["rate_continuous"]) if "rate_continuous" in overburden_history else [0.0]
        density_continuous_g_cm3 = np.atleast_1d(overburden_history["density_continuous_g_cm3"]) if "density_continuous_g_cm3" in overburden_history else [1.0]

        time_discrete = np.atleast_1d(overburden_history["time_discrete_kyr"]) if "time_discrete_kyr" in overburden_history else [0.0]
        overburden_discrete = np.atleast_1d(overburden_history["overburden_discrete"]) if "overburden_discrete" in overburden_history else [0.0]
        density_discrete_g_cm3 = np.atleast_1d(overburden_history["density_discrete_g_cm3"]) if "density_discrete_g_cm3" in overburden_history else [1.0]

        times = np.linspace(0, self.total_age_kyr, 10000)
        overburdens = np.zeros_like(times) + initial_depth * initial_density_g_cm3

        for start, end, rate, density in zip(start_time_continuous_kyr, end_time_continuous_kyr, rate_continuous, density_continuous_g_cm3):
            mask = (times >= start) & (times <= end)
            overburdens[mask] += density * rate * (times[mask] - start)
            mask_after_end = (times > end)
            overburdens[mask_after_end] += density * rate * (end - start)

        for t, o, d in zip(time_discrete, overburden_discrete, density_discrete_g_cm3):
            mask = (times >= t)
            overburdens[mask] += d * o

        if self.verbose>0:
            print("Interpolating overburden history...")

        return interp1d(times, overburdens, bounds_error=True)

    def _load_nuclear_data(self, filename):
        """
        Loads and caches nuclear data files (e.g., U238.dat, BindingEne.txt).

        Args:
            filename (str): The name of the nuclear data file to load.

        Returns:
            np.ndarray: The loaded data as a NumPy array.
        """
        if filename in self._nuclear_data_cache:
            return self._nuclear_data_cache[filename]
        
        filepath = os.path.join(self.data_path,"nuclear_data", filename)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Nuclear data file not found: {filepath}")

        if "BindingEne.txt" in filename:
            cols_to_use = (0, 1, 2)
        elif "U238.dat" in filename:
            cols_to_use = (1, 2, 3)
        else:
            raise ValueError(f"Unknown nuclear data file format: {filename}")
            
        self._nuclear_data_cache[filename] = np.loadtxt(filepath, usecols=cols_to_use, unpack=True)
        return self._nuclear_data_cache[filename]
    
    def _load_srim_data(self, ion_z):
        """
        Loads and caches processed SRIM data for a given ion, creating an energy-to-range function.
        It checks for a pre-processed file; if not found, it processes the raw SRIM data.

        Args:
            ion_z (int): The atomic number (Z) of the ion.

        Returns:
            tuple: A tuple containing (e_to_x_func, e, dee_dx, den_dx, x) where:
                   e_to_x_func (function): Interpolation function from Energy [keV] to Track Length [µm].
                   e (np.ndarray): Energy array [keV].
                   dee_dx (np.ndarray): Electronic stopping power [keV/µm].
                   den_dx (np.ndarray): Nuclear stopping power [keV/µm].
                   x (np.ndarray): Track length array [µm].
        """
        if ion_z in self._srim_cache:
            return self._srim_cache[ion_z]

        ion_symbol = element(ion_z).symbol
        raw_srim_filename = f"{element(ion_z).name} in {self.composition}.txt"
        raw_srim_filepath = os.path.join(self.data_path, "SRIM_data", self.name, raw_srim_filename)    
        processed_srim_dir = os.path.join(self.data_path, "SRIM_data", self.name)
        processed_srim_filepath = os.path.join(processed_srim_dir, f"{ion_symbol}-{self.shortname}.txt")

        e_kev, dee_dx, den_dx, length_um = (None, None, None, None)

        if os.path.exists(processed_srim_filepath):
            e_kev, dee_dx, den_dx, length_um = np.loadtxt(processed_srim_filepath, skiprows=1, unpack=True)
        else:
            if self.verbose>0:
                print(f"Pre-processed SRIM data for {ion_symbol} (Z={ion_z}) not found.")

            e_kev, dee_dx, den_dx, length_um = self._process_raw_srim_data(raw_srim_filepath)
            if e_kev is not None and length_um is not None:
                np.savetxt(processed_srim_filepath, np.column_stack((e_kev, dee_dx, den_dx, length_um)), header="Energy(keV)    dEe/dx(keV/micro_m)  dEn/dx(keV/micro_m)  x(micro_m)", fmt="%.6e")

        if e_kev is None or length_um is None:
            print(f"Warning: Could not load or process SRIM data for Z={ion_z}. Skipping.")
            return None, None, None, None, None
 
        unique_length, unique_indices = np.unique(length_um, return_index=True)
        unique_e = e_kev[unique_indices]

        if len(unique_e) < 2:
            print(f"Error: Not enough unique data points for interpolation for Z={ion_z}. Skipping.")
            return None, None, None, None, None

        e_to_x_func = interp1d(unique_e, unique_length, bounds_error=False, fill_value="extrapolate")
        self._srim_cache[ion_z] = (e_to_x_func, unique_e, dee_dx[unique_indices], den_dx[unique_indices], unique_length)
        return e_to_x_func, unique_e, dee_dx[unique_indices], den_dx[unique_indices], unique_length
    
    def _process_raw_srim_data(self, raw_srim_filepath):
        """
        Processes a raw SRIM data file, handling unit conversions.

        Args:
            raw_srim_filepath (str): The full path to the raw SRIM output file.

        Returns:
            tuple: A tuple containing (e_kev, dee_dx, den_dx, length_um) where:
                   e_kev (np.ndarray): Energy array [keV].
                   dee_dx (np.ndarray): Electronic stopping power [keV/µm].
                   den_dx (np.ndarray): Nuclear stopping power [keV/µm].
                   length_um (np.ndarray): Projected track length [µm].
        """
        if not os.path.exists(raw_srim_filepath):
            print(f"Error: Raw SRIM file not found at {raw_srim_filepath}.")
            return None, None, None, None

        def comma_to_dot(x):
            if isinstance(x, bytes):
                x = x.decode()
            return float(x.replace(',', '.'))

        data = np.genfromtxt(raw_srim_filepath, usecols=(0, 2, 3, 4, 6, 8), unpack=True, skip_header=(22+len(self.config['composition'].split('-'))), skip_footer=13, converters={0: comma_to_dot, 2: comma_to_dot, 3: comma_to_dot, 4: comma_to_dot, 6: comma_to_dot, 8: comma_to_dot})

        e_raw, dee_dx, den_dx, x_raw, y_raw, z_raw = data
        unit_e, unit_x = np.genfromtxt(raw_srim_filepath, dtype=str, skip_header=(22+len(self.config['composition'].split('-'))), skip_footer=13, usecols=(1, 5), unpack=True)

        e_kev = np.zeros_like(e_raw)
        for j, unit in enumerate(unit_e):
            if unit == "eV": e_kev[j] = e_raw[j] * 1e-3
            elif unit == "MeV": e_kev[j] = e_raw[j] * 1e3
            else: e_kev[j] = e_raw[j]

        x_um, y_um, z_um = np.zeros_like(x_raw), np.zeros_like(y_raw), np.zeros_like(z_raw)
        for j, unit in enumerate(unit_x):
            if unit == "A":
                x_um[j], y_um[j], z_um[j] = x_raw[j] * 1e-4, y_raw[j] * 1e-4, z_raw[j] * 1e-4
            else:
                x_um[j], y_um[j], z_um[j] = x_raw[j], y_raw[j] * 1e-4, z_raw[j] * 1e-4

        length_um = np.sqrt(x_um**2 + y_um**2 + z_um**2)

        return e_kev, dee_dx, den_dx, length_um
    
    def _radiogenic_spectrum(self):
        """
        Calculates the absolute neutron flux from Spontaneous fission.
        Based on Cranberg et al., Phys. Rev. 103, 662 (1956), https://www.wise-uranium.org/ranc.html

        Returns:
            interpolator (np.interp1d): Interpolator of energy vs fluxes in g^-1 s^-1 GeV^-1.
        """
        sf_yield_per_gu_s = 0.1353
        sf_yield_per_g_s = self.config['uranium_concentration_g_g'] * sf_yield_per_gu_s
        energies = np.logspace(-6, -1, 10000)

        a, b = 0.00065, 3700.0
        watt_shape = np.exp(-energies / a) * np.sinh(np.sqrt(b * energies))
        
        sf_flux = watt_shape * (sf_yield_per_g_s/ np.trapezoid(watt_shape, energies))

        interpolator = interp1d(energies, sf_flux, bounds_error=False, fill_value='extrapolate')

        return interpolator
    
    def _alpha_n_spectrum(self):
        """
        Calculates the absolute neutron flux from alpha,n reactions.
        Based on Kudryavtsev at al., SciPost Phys. Proc. 12, 2023 (SOURCES4)
        
        Returns:
            interpolator (np.interp1d): Interpolator of energy vs fluxes in g^-1 s^-1 GeV^-1.
        """
        an_yield_per_gu_s = self.config["an_coeff_n_nalpha"] * 8 * 1.245 * 1e-5 * 1e9 / 1e6

        an_yield_per_g_s = an_yield_per_gu_s * self.config['uranium_concentration_g_g']

        alpha, beta_gev, cutoff_e_gev = 2, 900, 0.0065

        energies_gev = np.logspace(-4, -0.5, 10000)

        an_shape_gev = (energies_gev**alpha)* np.exp(-beta_gev*energies_gev)
        an_shape_gev *= 1.0 / (1.0 + np.exp((energies_gev - cutoff_e_gev) / 0.0005))

        an_flux_gev = an_shape_gev * (an_yield_per_g_s / np.trapezoid(an_shape_gev, energies_gev))

        interpolator = interp1d(energies_gev, an_flux_gev, bounds_error=False, fill_value='extrapolate')

        return interpolator

    def _process_background_neutron_geant4_data(self, background_type, energy_bins_gev, total_simulated_particles=1e4):
        """
        Processes raw Geant4 data for neutrons from spontaneous fission or from alpha,n reactions, creating a normalized recoil spectrum file.

        Args:
            background_type (str): The specific background to compute for.
            energy_bins_gev (np.ndarray): The energy bin edges [GeV].
            total_simulated_particles (float, optional): Number of particles per Geant4 run. Defaults to 1e4.
        """
        
        if background_type == 'fission_n':
            flux_interpolator = self.radiogenic_spectrum
        elif background_type == 'alpha_n':
            flux_interpolator = self.alpha_n_spectrum

        all_fragments = self._get_all_fragments(energy_bins_gev[:-1], 'neutron')

        geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_neutron")
                
        all_recoil_spectra = {}

        fragment_spectra = {frag: np.zeros(len(RECOIL_ENERGY_BINS_MEV) - 1) for frag in all_fragments}

        for i in range(len(energy_bins_gev) - 1):
            e_min = energy_bins_gev[i]
            e_max = energy_bins_gev[i+1]

            N_E = 50
            e_vals = np.linspace(e_min, e_max, N_E)

            flux_grid = flux_interpolator(e_vals)

            weight = np.trapezoid(flux_grid, e_vals)

            filepath = os.path.join(geant4_input_dir, f"outNuclei_{e_min:.6f}.txt")
            if not os.path.exists(filepath): continue

            try:
                df = pd.read_csv(filepath, sep=r'\s+', header=None, usecols=[0, 2], 
                                names=['name', 'rec_e'], dtype={'name': str, 'rec_e': float})
                names = df['name'].values
                rec_energies = df['rec_e'].values
            except pd.errors.EmptyDataError:
                continue

            valid_mask = np.isin(names, list(fragment_spectra.keys()))
            valid_names = names[valid_mask]
            valid_rec = rec_energies[valid_mask]

            unique_frags = np.unique(valid_names)

            for frag in unique_frags:
                frag_mask = (valid_names == frag)
                frag_rec_energies = valid_rec[frag_mask]

                counts, _ = np.histogram(frag_rec_energies, bins=RECOIL_ENERGY_BINS_MEV)
                fragment_spectra[frag] += counts * weight

        all_recoil_spectra.update(fragment_spectra)

        bin_widths_mev = np.diff(RECOIL_ENERGY_BINS_MEV)
        norm_factor = (bin_widths_mev * total_simulated_particles * MYR_PER_SECOND * 1e-3)

        output_dir = os.path.join(self.data_path, "processed_recoils", self.name)
        os.makedirs(output_dir, exist_ok=True)
        output_filepath = os.path.join(output_dir, f"{background_type}.npz")

        normalized_spectra = {}
        for name, spectrum in all_recoil_spectra.items():
            normalized_spectra[name] = np.divide(spectrum, norm_factor, out=np.zeros_like(spectrum), where=norm_factor!=0)

        np.savez(output_filepath, Er_bins=RECOIL_ENERGY_BINS_MEV, **normalized_spectra)
        if self.verbose>1:
            print(f"    - Saved processed data to {output_filepath}")
        
    def calculate_background_neutron_spectrum(self, 
        x_bins, 
        energy_bins_gev, 
        background_types=['fission_n', 'alpha_n'],
        total_simulated_particles=1e4, 
        ): 
        """
        Calculates the neutron-induced recoil spectrum for a specific background type.

        Args:
            background_type (str): The specific background to compute for.
            energy_bins_gev (np.ndarray): The energy bin edges [GeV].
            total_simulated_particles (float, optional): Number of particles per Geant4 run. Defaults to 1e4.
        """

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0

        sum_drdx = np.zeros_like(x_mids)

        for background_type in background_types:
            filepath = os.path.join(self.data_path, "processed_recoils", self.name, f"{background_type}.npz")
            if not os.path.exists(filepath):
                self._process_background_neutron_geant4_data(background_type, energy_bins_gev, total_simulated_particles)

            recoil_data = np.load(filepath)

            drdx = self._convert_recoil_to_track_spectrum(x_bins, recoil_data=recoil_data, energy_bins_gev=energy_bins_gev, species='neutron')

            sum_drdx += drdx['total']

        return sum_drdx

    def integrate_background_neutron_spectrum(
        self, 
        x_bins, 
        energy_bins_gev,  
        sample_mass_kg,
        background_types=['fission_n', 'alpha_n'],
        total_simulated_particles=1e4, 
        x_grid=TRACK_LENGTH_BINS_NM, 
        ):
        """
        Integrates the neutron-induced recoil spectrum over specified background types to obtain total track counts.
        
        Args:
            x_bins (np.ndarray): The bin edges for the output track length spectrum [nm].
            energy_bins_gev (np.ndarray): The energy bin edges [GeV].
            total_exposure_kyr (float): The total exposure time [kyr].
            sample_mass_kg (float): The mass of the mineral sample [kg].
            background_types (list, optional): List of background types to include. Defaults to ['fission', 'alpha_n'].
            total_simulated_particles (float, optional): Number of particles per Geant4 run. Defaults to 1e4.
            x_grid (np.ndarray, optional): The grid for track length spectrum calculation. Defaults to TRACK_LENGTH_BINS_NM.
        Returns:
            np.ndarray: The total track counts in each bin.
        """
        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
        x_mids_grid = x_grid[:-1] + np.diff(x_grid) / 2.0

        sum_drdx = self.calculate_background_neutron_spectrum(
            x_grid, 
            energy_bins_gev, 
            background_types,
            total_simulated_particles, 
            )

        sum_drdx *= self.total_age_kyr * 1e-3 * sample_mass_kg

        total_tracks_interp  = interp1d(x_mids_grid, np.array(sum_drdx),  bounds_error=False, fill_value='extrapolate')

        total_tracks = np.array([quad(total_tracks_interp, x_bins[i], x_bins[i+1])[0] for i in range(len(x_mids))])

        return total_tracks

    def calculate_nu_spectrum(self, x_bins=TRACK_LENGTH_BINS_NM, flux_name='all'):
        """
        Calculates the differential track rate from neutrino sources via CEvNS.

        Args:
            x_bins (np.ndarray): The bin edges for the output track length spectrum [nm].
            flux_name (str, optional): The neutrino flux source as defined in WIMpy. Defaults to 'all'.

        Returns:
            np.ndarray: The differential track rate (dR/dx) [events/kg/Myr/nm].
        """
        if self.verbose>0:
            print(f"Calculating neutrino background for source: {flux_name}...")
        try:        
            from WIMpy import DMUtils as DMU
        except ImportError as e:
            raise RuntimeError(
                "The `WIMpy` package is required for this function. "
                "Please install to use neutrino spectrum calculations."
            ) from e
        
        x_mid = x_bins[:-1] + np.diff(x_bins) / 2.0
        dRdx = np.zeros_like(x_mid)

        for i, nuc_name in enumerate(self.config['nuclei']):
            if nuc_name != "H":
                ion_z = element(nuc_name).atomic_number
                srim_func, e, dee_dx, den_dx, x = self._load_srim_data(ion_z)
                if not srim_func: continue

                sorted_indices = np.argsort(x)
            
                x_to_e_func = interp1d(x[sorted_indices]*1e3, e[sorted_indices], bounds_error=False, fill_value=0.0)
                x_to_dedx_func = interp1d(x[sorted_indices]*1e3, (dee_dx[sorted_indices]+den_dx[sorted_indices])*1e-3, bounds_error=False, fill_value=0.0)

                dRdE_kev = np.vectorize(DMU.dRdE_CEvNS)(x_to_e_func(x_mid), element(nuc_name).protons, element(nuc_name).neutrons, flux_name=flux_name)

                dRdx += self.config['stoich'][i] * dRdE_kev * np.abs(x_to_dedx_func(x_mid))
        
        return dRdx * 365 * 1e6
    
    def integrate_nu_spectrum(self, x_bins, sample_mass, flux_name="all", x_grid=TRACK_LENGTH_BINS_NM):

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
        x_mids_grid = x_grid[:-1] + np.diff(x_grid) / 2.0

        drdx = self.calculate_nu_spectrum(x_grid, flux_name) * self.total_age_kyr * sample_mass * 1e-3

        total_tracks_interp  = interp1d(x_mids_grid, np.array(drdx),  bounds_error=False, fill_value='extrapolate')

        total_tracks = [quad(total_tracks_interp, x_bins[i], x_bins[i+1])[0] for i in range(len(x_mids))]

        return total_tracks
    
    def calculate_fission_spectrum(self, x_bins=TRACK_LENGTH_BINS_NM):
        """
        Calculates the differential track rate from U-238 spontaneous fission.

        Args:
            x_bins (np.ndarray): The bin edges for the output track length spectrum [nm].

        Returns:
            np.ndarray: The differential track rate (dR/dx) [events/kg/Myr/nm].
        """
        if self.verbose>0:
            print("Calculating spontaneous fission background...")
        
        Z_fission, A_fission, _ = self._load_nuclear_data("U238.dat")
        Z_bind, A_bind, B_bind = self._load_nuclear_data("BindingEne.txt")
        binding_map = {(int(z), int(a)): b for z, a, b in zip(Z_bind, A_bind, B_bind)}

        B0_U238 = 7.570126
        tau_U238_fission = 1.21e11  # Myr
        mass_U238 = 238.0
        u_concentration = self.config["uranium_concentration_g_g"]

        fission_rate_factor = (u_concentration * (AVOGADRO_NUMBER / mass_U238) * 1e3) / tau_U238_fission
        
        total_track_lengths_nm = []
        num_events = int(len(Z_fission) / 3)

        for i in range(num_events):
            z1, a1 = int(Z_fission[3*i + 1]), int(A_fission[3*i + 1])
            z2, a2 = int(Z_fission[3*i + 2]), int(A_fission[3*i + 2])

            b1 = binding_map.get((z1, a1), 0)
            b2 = binding_map.get((z2, a2), 0)
            if b1 == 0 or b2 == 0: continue

            M0 = 92 * PROTON_MASS_MEV + (238 - 92) * NEUTRON_MASS_MEV - B0_U238 * 238
            m1 = z1 * PROTON_MASS_MEV + (a1 - z1) * NEUTRON_MASS_MEV - b1 * a1
            m2 = z2 * PROTON_MASS_MEV + (a2 - z2) * NEUTRON_MASS_MEV - b2 * a2

            Ek1_MeV = (M0**2 + m1**2 - m2**2) / (2 * M0) - m1
            Ek2_MeV = (M0**2 + m2**2 - m1**2) / (2 * M0) - m2

            srim_func1, _, _, _, _ = self._load_srim_data(z1)
            srim_func2, _, _, _, _ = self._load_srim_data(z2)

            if srim_func1 and srim_func2:
                track1 = srim_func1(Ek1_MeV * 1e3) * 1e3
                track2 = srim_func2(Ek2_MeV * 1e3) * 1e3
                total_track_lengths_nm.append(track1 + track2)
                
        counts, bin_edges = np.histogram(total_track_lengths_nm, bins=x_bins)
        bin_widths = np.diff(bin_edges)
        
        dRdx = (counts / num_events) * fission_rate_factor / bin_widths
        
        return dRdx
    
    def integrate_fission_spectrum(self, x_bins, sample_mass, x_grid=TRACK_LENGTH_BINS_NM):

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
        x_mids_grid = x_grid[:-1] + np.diff(x_grid) / 2.0

        drdx = self.calculate_fission_spectrum(x_grid) * self.total_age_kyr * sample_mass * 1e-3

        total_tracks_interp  = interp1d(x_mids_grid, np.array(drdx),  bounds_error=False, fill_value='extrapolate')

        total_tracks = np.asarray([quad(total_tracks_interp, x_bins[i], x_bins[i+1])[0] for i in range(len(x_mids))])

        return total_tracks

    def _load_depth_interpolators(self, species='mu-'):
        """
        Loads depth interpolators for maximum penetration depth and mean width
        from Geant4 simulation data for a given particle species.

        Args:
            species (str, optional): The particle species to simulate ('mu+', 'mu-', or 'neutron'). Defaults to 'mu-'.
        """
        pen = []
        width = []
        attenuation = []

        pri_energies = np.logspace(-3, 3, 10)

        tab_species = 'mu-' if species == 'mu-' or species == 'mu+' else 'neutron'

        cylinder = 500. if species == 'mu-' or species == 'mu+' else 100.

        for energy_name in pri_energies:

            filepath = os.path.join(self.data_path, "Geant4_data", f"StdRock_{tab_species}", f"outNuclei_{energy_name:.6f}.txt")
            if not os.path.exists(filepath):
                print(f"Warning: Geant4 data file not found: {filepath}")
                continue

            energies, depth, rem_energy = np.loadtxt(filepath, usecols=(2, 3, 5), dtype = str, unpack=True)
            energies = energies.astype(float)
            depth = cylinder - depth.astype(float)/1000.
            rem_energy = rem_energy.astype(float)

            if tab_species == 'mu-':
                pen.append(depth[rem_energy == 0].mean()*2.65)
                width.append(depth[rem_energy == 0].std()*2.65)
                
            else:
                steps = np.linspace(0, 1, 50)
                midsteps = steps[:-1] + np.diff(steps) / 2
                counts = np.histogram(depth*2.65, bins = steps)[0]
                slope = np.polyfit(midsteps, np.log(counts), 1)[0]
                attenuation.append(-1./slope)
        
        self._depth_interpolators[species] = {}

        if tab_species == 'mu-':

            self._maxdepth_interp = interp1d(pri_energies[:-1], pen[:-1], kind='linear', fill_value='extrapolate', bounds_error=False)
            self._meanwidth_interp = interp1d(pri_energies[:-1], width[:-1], kind='linear', fill_value='extrapolate', bounds_error=False)

            self._depth_interpolators[species]['maxdepth'] = self._clipped_maxdepth
            self._depth_interpolators[species]['meanwidth'] = self._clipped_meanwidth

        else:
            self._depth_interpolators[species]['attenuation'] = interp1d(pri_energies, attenuation, kind='linear', fill_value='extrapolate', bounds_error=False)
    
    def _clipped_maxdepth(self, x):
        """Clipped maximum penetration depth (pickleable instance method)."""
        return np.clip(self._maxdepth_interp(x), 0.5e-3, np.inf)

    def _clipped_meanwidth(self, x):
        """Clipped mean width (pickleable instance method)."""
        return np.clip(self._meanwidth_interp(x), 1.e-4, np.inf)

    def _flux_time_kyr(self, t_kyr):
        """
        Converts this class's local exposure-clock time (0 = start of data
        taking / formation, total_age_kyr = present) into the absolute
        timeline used by FluxHistory (0 = present, negative = past). This
        is the single place that mapping happens; every flux lookup below
        goes through it.

        Args:
            t_kyr (float): Local time [kyr], 0 (start of data taking) to
                total_age_kyr (present).
        Returns:
            float: The corresponding time [kyr] on the FluxHistory timeline.
        """
        return t_kyr - self.total_age_kyr
    
    def _get_local_neutron_flux_batch(
        self, target_depth_array, t_kyr_array, energy_bins_gev,
        total_simulated_particles=1e4, species_list=('mu-', 'mu+', 'neutron'),
        n_depth_bins=50,
    ):
        """
        Args:
            target_depth_array (np.ndarray): Physical depth [cm] to which each
                timestep's overburden corresponds (depth_mwe / density), shape (T,).
            t_kyr_array (np.ndarray): Local exposure-clock times [kyr], shape (T,).
            n_depth_bins (int): Resolution of the SHARED depth grid (matches the
                original hardcoded 50).

        Returns:
            tuple:
                slice_yield (np.ndarray): shape (T, n_depth_bins, len(energy_bins_gev))
                depth_bins (np.ndarray): the shared physical-depth bin edges [cm],
                    shape (n_depth_bins + 1,) -- sized to max(target_depth_array).
        """
        t_kyr_array = np.atleast_1d(t_kyr_array)
        target_depth_array = np.atleast_1d(target_depth_array)
        n_t = len(t_kyr_array)

        max_depth = np.max(target_depth_array)
        depth_bins = np.linspace(0, max_depth, n_depth_bins + 1)

        internal_edges = 0.5 * (energy_bins_gev[:-1] + energy_bins_gev[1:])
        first_edge = energy_bins_gev[0] - (internal_edges[0] - energy_bins_gev[0])
        last_edge = energy_bins_gev[-1] + (energy_bins_gev[-1] - internal_edges[-1])
        edges = np.concatenate(([first_edge], internal_edges, [last_edge]))

        t_kyr_flux_array = self._flux_time_kyr(t_kyr_array)  # (T,)

        slice_yield = np.zeros((n_t, n_depth_bins, len(energy_bins_gev)))

        for species in species_list:
            geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_{species}")

            initial_flux = np.stack(
                [self.flux_history.flux(species, energy_bins_gev, t) for t in t_kyr_flux_array],
                axis=0,
            )
            accumulated_yield = np.zeros((n_t, len(edges) - 1))

            counts = np.zeros((len(energy_bins_gev), n_depth_bins, len(edges) - 1))

            for i, e in enumerate(energy_bins_gev):
                filepath = os.path.join(geant4_input_dir, f"outNuclei_{e:.6f}.txt")
                if not os.path.exists(filepath):
                    continue

                df = pd.read_csv(
                    filepath, sep=r'\s+', header=None, usecols=[0, 2, 3],
                    names=['name', 'rec_e', 'depth'],
                    dtype={'name': 'category', 'rec_e': 'float32', 'depth': 'float32'},
                    engine='c',
                )
                df_n = df[df['name'] == 'neutron']

                converted_depths = 100.0 - (df_n['depth'].values * 1e-3)
                rec_energies_gev = df_n['rec_e'].values * 1e-3
                H, _, _ = np.histogram2d(converted_depths, rec_energies_gev, bins=[depth_bins, edges])
                counts[i, :, :] += H / total_simulated_particles

            for k in range(n_depth_bins):
                current_flux = initial_flux + accumulated_yield          # (T, n_energy)
                produced_at_k = current_flux @ counts[:, k, :]            # (T, n_energy)@(n_energy, n_edges-1)
                accumulated_yield += produced_at_k
                slice_yield[:, k, :] += produced_at_k

        return slice_yield, depth_bins

    def _get_all_fragments(self, energy_names_gev, species='mu-'):
        """
        Dynamically determines the list of all nuclear fragments from Geant4 output files.

        Args:
            energy_names_gev (list): List of energies used in Geant4 filenames.
            species (str, optional): The particle species to simulate ('mu+', 'mu-', or 'neutron'). Defaults to 'mu-'.


        Returns:
            list: A sorted list of unique fragment symbols.
        """
        geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_{species}")
        all_fragments = set()
        
        for energy_name in energy_names_gev:
            filepath = os.path.join(geant4_input_dir, f"outNuclei_{energy_name:.6f}.txt")
            if not os.path.exists(filepath): continue
            elif os.path.getsize(filepath) > 0:
                names = np.loadtxt(filepath, usecols=0, dtype=str, ndmin=1)

                for name in names:
                    if name not in ['He3', 'He4', 'He5', 'He6', 'He7', 'He8', 'alpha', 'proton']:
                        all_fragments.add(name)
        return sorted(list(all_fragments))

    def _process_geant4_data(
        self,
        t_kyr_array,
        energy_bins_gev,
        total_simulated_particles=1e4,
        target_thickness_mm=TYPICAL_DEPTH_MM,
        species='mu-',
    ):
        """
        Vectorized version: processes raw Geant4 data for an ARRAY of times at
        once, reading/histogramming each energy-bin file exactly once and
        reusing it across all requested timesteps.

        Args:
            t_kyr_array (np.ndarray): Local exposure-clock times [kyr], shape (T,).
            energy_bins_gev (np.ndarray): Energy bin edges [GeV].
            total_simulated_particles (float): Particles per Geant4 run.
            target_thickness_mm (float): Target slice thickness [mm].
            species (str): 'mu-', 'mu+', or 'neutron'.

        Returns:
            dict: {
                't_kyr': t_kyr_array,
                'depth_mwe': depth_mwe_array,           # (T,)
                'Er_bins': RECOIL_ENERGY_BINS_MEV,
                <fragment_name>: array of shape (T, n_recoil_bins), ...
            }
            Also writes this dict to a single cache file for the whole batch.
        """
        if self.flux_history is None:
            raise ValueError(
                "flux_history not initialized. Call integrate_particle_signal_spectrum_parallel "
                "first, or set self.flux_history directly."
            )
        if not self._depth_interpolators.get(species):
            raise ValueError(f"Depth interpolators not initialized for species {species}.")

        t_kyr_array = np.atleast_1d(np.asarray(t_kyr_array, dtype=float))
        n_t = t_kyr_array.shape[0]

        depth_mwe_array = self._overburden_interpolator(t_kyr_array)
        t_kyr_flux_array = self._flux_time_kyr(t_kyr_array)

        if species in ('mu-', 'mu+'):
            maxdepth = self._depth_interpolators[species]['maxdepth']
            meanwidth = self._depth_interpolators[species]['meanwidth']
        else:
            slope = self._depth_interpolators[species]['attenuation']

        all_fragments = self._get_all_fragments(energy_bins_gev[:-1], species)
        geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_{species}")

        n_recoil_bins = len(RECOIL_ENERGY_BINS_MEV) - 1
        fragment_spectra = {frag: np.zeros((n_t, n_recoil_bins)) for frag in all_fragments}

        target_thickness_mwe = target_thickness_mm * 0.001 * self.config['density_g_cm3']

        for i in range(len(energy_bins_gev) - 1):
            e_min = energy_bins_gev[i]
            e_max = energy_bins_gev[i + 1]

            N_E, N_C = 50, 50
            e_vals = np.linspace(e_min, e_max, N_E)
            c_vals = np.linspace(0.01, 1.0, N_C)
            E_grid, C_grid = np.meshgrid(e_vals, c_vals)

            z_min_grid = depth_mwe_array[:, None, None] / C_grid[None, :, :]
            z_max_grid = z_min_grid + target_thickness_mwe / C_grid[None, :, :]

            flux_grid = self.flux_history.get_map(species, t_kyr_flux_array, E_grid)[2].reshape(t_kyr_flux_array.size, *E_grid.shape)

            if species in ('mu-', 'mu+'):
                D_grid = maxdepth(E_grid)
                W_grid = meanwidth(E_grid) 

                effective_min = np.maximum(z_min_grid, D_grid)
                top_A = np.minimum(z_max_grid, D_grid)
                val_A = np.maximum(0.0, top_A - z_min_grid)

                mask = z_max_grid > effective_min
                val_B = np.where(
                    mask,
                    W_grid * (
                        np.exp(-(np.maximum(effective_min, D_grid) - D_grid) / W_grid)
                        - np.exp(-(z_max_grid - D_grid) / W_grid)
                    ),
                    0.0,
                )

                prob_tail_grid = (val_A + val_B) / (D_grid + W_grid)
                integrand_tail = flux_grid * prob_tail_grid

                weight_elastic_t = np.trapezoid(
                    np.trapezoid(integrand_tail, e_vals, axis=-1), c_vals, axis=-1
                )

                if species == 'mu-':
                    prob_peak_grid = 0.5 * (
                        erf((z_max_grid - D_grid) / (np.sqrt(2) * W_grid))
                        - erf((z_min_grid - D_grid) / (np.sqrt(2) * W_grid))
                    )
                    integrand_peak = flux_grid * prob_peak_grid
                    weight_peak_t = np.trapezoid(
                        np.trapezoid(integrand_peak, e_vals, axis=-1), c_vals, axis=-1
                    )
            else:
                lambda_grid = slope(E_grid)
                prob_attenuation_grid = (
                    np.exp(-z_min_grid / lambda_grid) - np.exp(-z_max_grid / lambda_grid)
                )
                integrand = flux_grid * prob_attenuation_grid
                weight_elastic_t = np.trapezoid(
                    np.trapezoid(integrand, e_vals, axis=-1), c_vals, axis=-1
                )

            filepath = os.path.join(geant4_input_dir, f"outNuclei_{e_min:.6f}.txt")
            if not os.path.exists(filepath):
                continue
            try:
                df = pd.read_csv(
                    filepath, sep=r'\s+', header=None, usecols=[0, 2, 5],
                    names=['name', 'rec_e', 'rem_e'],
                    dtype={'name': str, 'rec_e': float, 'rem_e': float},
                )
                names = df['name'].values
                rec_energies = df['rec_e'].values
                rem_energies = df['rem_e'].values
            except pd.errors.EmptyDataError:
                continue

            valid_mask = np.isin(names, list(fragment_spectra.keys()))
            valid_names = names[valid_mask]
            valid_rec = rec_energies[valid_mask]
            valid_rem = rem_energies[valid_mask]
            unique_frags = np.unique(valid_names)

            for frag in unique_frags:
                frag_mask = (valid_names == frag)
                frag_rec_energies = valid_rec[frag_mask]
                frag_rem_energies = valid_rem[frag_mask]

                if species == 'mu-':
                    mask_peak = (frag_rem_energies == 0.0)
                    rec_peak = frag_rec_energies[mask_peak]
                    rec_elastic = frag_rec_energies[~mask_peak]

                    if len(rec_peak) > 0:
                        counts_peak, _ = np.histogram(rec_peak, bins=RECOIL_ENERGY_BINS_MEV)
                        fragment_spectra[frag] += np.outer(weight_peak_t, counts_peak)

                    if len(rec_elastic) > 0:
                        counts_elastic, _ = np.histogram(rec_elastic, bins=RECOIL_ENERGY_BINS_MEV)
                        fragment_spectra[frag] += np.outer(weight_elastic_t, counts_elastic)
                else:
                    counts, _ = np.histogram(frag_rec_energies, bins=RECOIL_ENERGY_BINS_MEV)
                    fragment_spectra[frag] += np.outer(weight_elastic_t, counts)

        bin_widths_mev = np.diff(RECOIL_ENERGY_BINS_MEV)
        norm_factor = (
            bin_widths_mev * target_thickness_mm * self.config['density_g_cm3']
            * total_simulated_particles * MYR_PER_SECOND
        )

        normalized_spectra = {}
        for name, spectrum in fragment_spectra.items():
            if name == 'neutron':
                denom = bin_widths_mev * total_simulated_particles
                normalized_spectra[name] = np.divide(
                    spectrum, denom[None, :], out=np.zeros_like(spectrum), where=denom[None, :] != 0
                )
            else:
                normalized_spectra[name] = np.divide(
                    spectrum, norm_factor[None, :], out=np.zeros_like(spectrum),
                    where=norm_factor[None, :] != 0,
                )

        scenario_name = self.flux_history.baseline.name + "_" + "_".join(
            event["template"].name + "_" + str(-event["start_time_kyr"]) + "kyr"
            for event in self.flux_history.events
        )

        output_dir = os.path.join(self.data_path, "processed_recoils", self.name, species)
        os.makedirs(output_dir, exist_ok=True)
        output_filepath = os.path.join(output_dir, f"{scenario_name}_depth{depth_mwe_array[0]:.1f}_{depth_mwe_array[-1]:.1f}mwe.npz")

        np.savez(
            output_filepath,
            t_kyr=t_kyr_array,
            depth_mwe=depth_mwe_array,
            Er_bins=RECOIL_ENERGY_BINS_MEV,
            **normalized_spectra,
        )
        if self.verbose > 1:
            print(f"    - Saved batched processed data ({n_t} timesteps) to {output_filepath}")

        return {
            't_kyr': t_kyr_array,
            'depth_mwe': depth_mwe_array,
            'Er_bins': RECOIL_ENERGY_BINS_MEV,
            **normalized_spectra,
        }


    def _process_secondary_geant4_data(
        self, t_kyr_array, energy_bins_gev,
        total_simulated_particles=1e4, target_thickness_mm=TYPICAL_DEPTH_MM,
        secondary_neutrons_species=('mu-', 'mu+', 'neutron'),
    ):
        """
        Vectorized replacement: processes secondary-neutron Geant4 data for an
        ARRAY of times at once. See _get_local_neutron_flux_batch docstring for
        how the time-dependent depth range (via overburden) is handled.
        """
        if not self._depth_interpolators.get('neutron'):
            raise ValueError("Depth interpolators not initialized for neutrons.")

        slope = self._depth_interpolators['neutron']['attenuation']
        all_fragments = self._get_all_fragments(energy_bins_gev[:-1], species='neutron')
        geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_neutron")

        t_kyr_array = np.atleast_1d(np.asarray(t_kyr_array, dtype=float))
        n_t = len(t_kyr_array)

        density = self.config['density_g_cm3']
        depth_mwe_array = self._overburden_interpolator(t_kyr_array)
        target_depth_array = depth_mwe_array / density

        slice_yield, depth_bins = self._get_local_neutron_flux_batch(
            target_depth_array, t_kyr_array, energy_bins_gev,
            total_simulated_particles, secondary_neutrons_species,
        )

        eff_depth_bins = depth_bins * density
        eff_depth_mids = eff_depth_bins[:-1] + np.diff(eff_depth_bins) / 2.0
        depth_bin_widths_mwe = np.diff(eff_depth_bins)
        n_depth_bins = len(eff_depth_mids)

        n_recoil_bins = len(RECOIL_ENERGY_BINS_MEV) - 1
        fragment_spectra = {frag: np.zeros((n_t, n_recoil_bins)) for frag in all_fragments}
        target_thickness_mwe = target_thickness_mm * 0.001 * density

        for i in range(len(energy_bins_gev) - 1):
            e_min = energy_bins_gev[i]
            e_max = energy_bins_gev[i + 1]

            N_E, N_C = 50, 50
            e_vals = np.linspace(e_min, e_max, N_E)
            c_vals = np.linspace(0.01, 1.0, N_C)
            E_grid, C_grid = np.meshgrid(e_vals, c_vals)
            lambda_grid = slope(E_grid)

            weight_elastic_t = np.zeros(n_t)

            for j in range(n_depth_bins):
                eff_depth_t = depth_mwe_array - eff_depth_mids[j]   # (T,)
                valid_t = eff_depth_t > 0
                if not np.any(valid_t):
                    continue

                z_min_grid = eff_depth_t[valid_t][:, None, None] / C_grid[None, :, :]
                z_max_grid = z_min_grid + target_thickness_mwe / C_grid[None, :, :]
                prob_attenuation_grid = (
                    np.exp(-z_min_grid / lambda_grid) - np.exp(-z_max_grid / lambda_grid)
                )
                trapz_t = np.trapezoid(
                    np.trapezoid(prob_attenuation_grid, e_vals, axis=-1), c_vals, axis=-1
                )

                contribution = np.zeros(n_t)
                contribution[valid_t] = depth_bin_widths_mwe[j] * slice_yield[valid_t, j, i] * trapz_t
                weight_elastic_t += contribution

            filepath = os.path.join(geant4_input_dir, f"outNuclei_{e_min:.6f}.txt")
            if not os.path.exists(filepath):
                continue
            try:
                df = pd.read_csv(
                    filepath, sep=r'\s+', header=None, usecols=[0, 2],
                    names=['name', 'rec_e'], dtype={'name': str, 'rec_e': float},
                )
                names = df['name'].values
                rec_energies = df['rec_e'].values
            except pd.errors.EmptyDataError:
                continue

            valid_mask = np.isin(names, list(fragment_spectra.keys()))
            valid_names = names[valid_mask]
            valid_rec = rec_energies[valid_mask]

            for frag in np.unique(valid_names):
                frag_mask = (valid_names == frag)
                counts, _ = np.histogram(valid_rec[frag_mask], bins=RECOIL_ENERGY_BINS_MEV)
                fragment_spectra[frag] += np.outer(weight_elastic_t, counts)

        bin_widths_mev = np.diff(RECOIL_ENERGY_BINS_MEV)
        norm_factor = (
            bin_widths_mev * target_thickness_mm * density
            * total_simulated_particles * MYR_PER_SECOND
        )

        normalized_spectra = {}
        for name, spectrum in fragment_spectra.items():
            if name == 'neutron':
                denom = bin_widths_mev * total_simulated_particles
                normalized_spectra[name] = np.divide(
                    spectrum, denom[None, :], out=np.zeros_like(spectrum), where=denom[None, :] != 0
                )
            else:
                normalized_spectra[name] = np.divide(
                    spectrum, norm_factor[None, :], out=np.zeros_like(spectrum),
                    where=norm_factor[None, :] != 0,
                )

        scenario_name = self.flux_history.baseline.name + "_" + "_".join(
            event["template"].name + "_" + str(-event["start_time_kyr"]) + "kyr"
            for event in self.flux_history.events
        )
        output_dir = os.path.join(self.data_path, "processed_recoils", self.name, 'secondary_neutron')
        os.makedirs(output_dir, exist_ok=True)
        output_filepath = os.path.join(output_dir, f"{scenario_name}_depth{depth_mwe_array[0]:.1f}_{depth_mwe_array[-1]:.1f}mwe.npz")

        np.savez(
            output_filepath,
            t_kyr=t_kyr_array,
            depth_mwe=depth_mwe_array,
            Er_bins=RECOIL_ENERGY_BINS_MEV,
            **normalized_spectra,
        )
        if self.verbose > 1:
            print(f"    - Saved batched secondary-neutron data ({n_t} timesteps) to {output_filepath}")

        return {
            't_kyr': t_kyr_array,
            'depth_mwe': depth_mwe_array,
            'Er_bins': RECOIL_ENERGY_BINS_MEV,
            **normalized_spectra,
        }


    def _convert_recoil_to_track_spectrum(self, x_bins, recoil_data, energy_bins_gev, species='mu-'):
        """
        Same as before, but recoil_data[fragment] may now have shape (T, n_recoil_bins)
        instead of (n_recoil_bins,). scipy's interp1d(..., axis=-1) interpolates each
        time-row independently in one call.
        """
        er_bins = recoil_data['Er_bins']
        er_mid_mev = er_bins[:-1] + np.diff(er_bins) / 2.0

        tab_species = 'neutron' if species == 'secondary_neutron' else species
        all_fragments = self._get_all_fragments(energy_bins_gev[:-1], tab_species)

        x_mid_nm = x_bins[:-1] + np.diff(x_bins) / 2.0

        sample_frag = next((f for f in all_fragments if f in recoil_data and f != 'neutron'), None)

        input_was_1d = sample_frag is not None and np.asarray(recoil_data[sample_frag]).ndim == 1

        n_t = recoil_data[sample_frag].shape[0] if sample_frag is not None and recoil_data[sample_frag].ndim == 2 else 1

        dRdx_by_nucleus = {}
        dRdx_total = np.zeros((n_t, len(x_bins) - 1))

        for nuclide_name in all_fragments:
            if nuclide_name not in recoil_data or nuclide_name == 'neutron':
                continue

            dRdEr_mev = np.atleast_2d(recoil_data[nuclide_name])  # (T, n_recoil_bins)

            dRdEr_interp = interp1d(
                er_mid_mev, dRdEr_mev, axis=-1, bounds_error=False, fill_value=0.0
            )

            nucleus_name = ''.join(filter(str.isalpha, nuclide_name))
            ion_z = element(nucleus_name).atomic_number
            srim_func, e, dee_dx, den_dx, x = self._load_srim_data(ion_z)
            if not srim_func:
                continue

            sorted_indices = np.argsort(x)
            x_to_e_func = interp1d(x[sorted_indices] * 1e3, e[sorted_indices] * 1e-3,
                                    bounds_error=False, fill_value=0.0)
            x_to_dedx_func = interp1d(x[sorted_indices] * 1e3,
                                    (dee_dx[sorted_indices] + den_dx[sorted_indices]) * 1e-6,
                                    bounds_error=False, fill_value=0.0)

            e_at_x = x_to_e_func(x_mid_nm)

            dRdx_nucleus = dRdEr_interp(e_at_x) * x_to_dedx_func(x_mid_nm)[None, :]
            dRdx_by_nucleus[nuclide_name] = dRdx_nucleus
            dRdx_total += dRdx_nucleus

        dRdx_by_nucleus["total"] = dRdx_total

        if input_was_1d:
            dRdx_by_nucleus = {name: arr[0] for name, arr in dRdx_by_nucleus.items()}
 

        return dRdx_by_nucleus


    def calculate_particle_signal_spectrum(
        self, x_bins, t_kyr_array, energy_bins_gev,
        total_simulated_particles=1e4, target_thickness_mm=TYPICAL_DEPTH_MM,
        species='mu-', nucleus="total",
    ):
        """
        Batch replacement for calculate_particle_signal_spectrum: takes the whole
        t_kyr_array at once, returns dR/dx with shape (T, n_x_bins) [or dict of
        those, per-nucleus].
        """

        depth_mwe_array = self._overburden_interpolator(t_kyr_array)

        scenario_name = self.flux_history.baseline.name + "_" + "_".join(
            event["template"].name + "_" + str(-event["start_time_kyr"]) + "kyr"
            for event in self.flux_history.events
        )
        filepath = os.path.join(
            self.data_path, "processed_recoils", self.name, species, f"{scenario_name}_depth{depth_mwe_array[0]:.1f}_{depth_mwe_array[-1]:.1f}mwe.npz"
        )

        recoil_data = None
        if os.path.exists(filepath):
            cached = np.load(filepath)
            if cached['t_kyr'].shape == t_kyr_array.shape and np.allclose(cached['t_kyr'], t_kyr_array):
                recoil_data = cached

        if recoil_data is None:
            if species == 'secondary_neutron':
                recoil_data = self._process_secondary_geant4_data(
                    t_kyr_array, energy_bins_gev, total_simulated_particles, target_thickness_mm,
                )
            else:
                recoil_data = self._process_geant4_data(
                    t_kyr_array, energy_bins_gev, total_simulated_particles, target_thickness_mm, species,
                )

        dRdx_at_depth = self._convert_recoil_to_track_spectrum(x_bins, recoil_data, energy_bins_gev, species)

        if nucleus == "total":
            return dRdx_at_depth["total"]          # (T, n_x_bins)
        elif nucleus == "all":
            return dRdx_at_depth
        else:
            return dRdx_at_depth[nucleus]


    def integrate_particle_signal_spectrum(
        self, x_bins, energy_bins_gev, sample_mass_kg,
        exposure_window_kyr=None, flux_history=None, overburden_history=None,
        steps=None, total_simulated_particles=1e4, target_thickness_mm=TYPICAL_DEPTH_MM,
        x_grid=TRACK_LENGTH_BINS_NM, species='mu-',
    ):
        """
        Same public signature/behavior as before, but internally does ONE
        vectorized batch call instead of a Pool over timesteps. Multiprocessing
        now happens one level up, over species (see integrate_all_particles).
        """
        if flux_history is not None:
            self.flux_history = flux_history
        elif self.flux_history is None:
            self.flux_history = FluxHistory(
                baseline={"kind": "Baseline"}, template_dir=os.path.join(self.data_path, "flux_data")
            )

        self._load_depth_interpolators('neutron' if species == 'secondary_neutron' else species)

        if overburden_history is not None:
            self._overburden_interpolator = self._interpolate_overburden_history(overburden_history)

        if exposure_window_kyr is None:
            exposure_window_kyr = self.total_age_kyr
        if isinstance(exposure_window_kyr, (int, float)):
            exposure_window_kyr = [0, exposure_window_kyr]

        if not steps:
            steps = len(self.flux_history.events) + int(
                (exposure_window_kyr[1] - exposure_window_kyr[0])
            )

        t_kyr_array = np.linspace(exposure_window_kyr[0], exposure_window_kyr[1], steps)

        dRdx_array = self.calculate_particle_signal_spectrum(
            x_grid, t_kyr_array, energy_bins_gev,
            total_simulated_particles, target_thickness_mm, species,
        )

        total_drdx = np.trapezoid(dRdx_array, t_kyr_array, axis=0)

        spectrum_density = total_drdx * sample_mass_kg * 1e-3

        internal_bin_widths = np.diff(x_grid)
        cumulative_counts = np.concatenate(([0], np.cumsum(spectrum_density * internal_bin_widths)))
        cdf_interp = interp1d(x_grid, cumulative_counts, kind='linear', bounds_error=False,
                            fill_value=(0, cumulative_counts[-1]))
        total_tracks = cdf_interp(x_bins[1:]) - cdf_interp(x_bins[:-1])

        return total_tracks


    def integrate_all_particles(
        self, x_bins, energy_bins_gev, sample_mass_kg,
        exposure_window_kyr=None, flux_history=None, overburden_history=None,
        steps=None, total_simulated_particles=1e4, target_thickness_mm=TYPICAL_DEPTH_MM,
        species_list=('mu-', 'mu+', 'neutron', 'secondary_neutron'),
    ):
        """
        Same public behavior, but the multiprocessing.Pool now parallelizes
        across species (each worker does one fully-vectorized-over-time run)
        instead of across timesteps.
        """

        shared_kwargs = dict(
            x_bins=x_bins, energy_bins_gev=energy_bins_gev, sample_mass_kg=sample_mass_kg,
            exposure_window_kyr=exposure_window_kyr, flux_history=flux_history,
            overburden_history=overburden_history, steps=steps,
            total_simulated_particles=total_simulated_particles,
            target_thickness_mm=target_thickness_mm,
        )

        tasks = [(self, dict(shared_kwargs, species=s)) for s in species_list]

        with Pool(processes=min(len(species_list), os.cpu_count() or 1)) as pool:
            results = list(tqdm(pool.imap(_species_worker, tasks), total=len(tasks)))

        total_tracks_by_species = dict(results)
        total_tracks_by_species['total'] = sum(total_tracks_by_species.values())

        return total_tracks_by_species


def _species_worker(self_and_args):
    self, kwargs = self_and_args
    species = kwargs.pop('species')
    if self.verbose > 1:
        print(f"Processing species: {species}")
    return species, self.integrate_particle_signal_spectrum(species=species, **kwargs)