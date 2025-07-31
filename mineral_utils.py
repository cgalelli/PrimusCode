# mineral_utils.py
#
# A refactored and streamlined utility module for paleo-detector analysis.
# This module contains the core Mineral class and all necessary functions
# for calculating signal and background track spectra.

import numpy as np
import os
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.integrate import quad, cumtrapz
from mendeleev import element

# --- Physical Constants ---
PROTON_MASS_MEV = 938.3
NEUTRON_MASS_MEV = 939.6
AVOGADRO_NUMBER = 6.022e23
SECONDS_PER_MYR = 60 * 60 * 24 * 365 * 1e6

# --- Binning setup ---

RECOIL_N_BINS= 101
RECOIL_ER_MIN_LOG_MEV= -2 
RECOIL_ER_MAX_LOG_MEV= 3
RECOIL_ENERGY_BINS_MEV = np.logspace(RECOIL_ER_MIN_LOG_MEV, RECOIL_ER_MAX_LOG_MEV, RECOIL_N_BINS)


# --- Utility Functions ---

def log_interp1d(xx, yy, kind='linear'):
    """Helper function for log-log interpolation."""
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interp1d(logx, logy, kind=kind, fill_value='extrapolate')
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

# --- Main Mineral Class ---

class Mineral:
    """
    A class to handle all physics calculations for a specific mineral.
    
    This class is initialized with a configuration dictionary and provides
    methods to calculate various signal and background components for
    paleo-detector analysis.
    """
    
    def __init__(self, mineral_config, data_path_prefix="Data/"):
        """
        Initializes the Mineral object.

        Args:
            mineral_config (dict): A dictionary containing all properties of the mineral.
            data_path_prefix (str): The relative path to the data directory.
        """
        self.config = mineral_config
        self.name = mineral_config['name']
        self.shortname = mineral_config["shortname"]
        self.composition = mineral_config['composition']
        self.data_path = data_path_prefix
        
        # --- Caches to store loaded data for efficiency ---
        self._srim_cache = {}
        self._recoil_cache = {}
        self._nuclear_data_cache = {}
        
        print(f"Initialized Mineral: {self.name}")

    def _load_nuclear_data(self, filename):
        """Loads and caches nuclear data files (e.g., U238.dat)."""
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
        Loads and caches processed SRIM data for a given ion.
        This function creates an interpolation function from Energy [keV] to Track Length [nm].
        It first checks for a pre-processed file; if not found, it processes the raw SRIM data
        and saves the result for future use.
        """
        if ion_z in self._srim_cache:
            return self._srim_cache[ion_z]

        # Determine the atomic symbol for naming processed files
        ion_symbol = element(ion_z).symbol

        # Define paths for raw and processed SRIM data
        raw_srim_filename = f"{ion_z} in {self.composition}.txt"
        raw_srim_filepath = os.path.join(self.data_path, "SRIM_data", self.name, raw_srim_filename)    
        processed_srim_dir = os.path.join(self.data_path, "SRIM_data", self.name)
        processed_srim_filepath = os.path.join(processed_srim_dir, f"{ion_symbol}-{self.shortname}.txt")

        e_kev = None
        length_um = None

        # Try to load pre-processed data first
        if os.path.exists(processed_srim_filepath):
            print(f"Loading pre-processed SRIM data for {ion_symbol} (Z={ion_z}) from {processed_srim_filepath}")
            e_kev, dee_dx, den_dx, length_um = np.loadtxt(processed_srim_filepath, skiprows=1, unpack=True)
        else:
            # If processed file doesn't exist, process the raw file
            print(f"Pre-processed SRIM data for {ion_symbol} (Z={ion_z}) not found. Processing raw file.")
            e_kev, dee_dx, den_dx, length_um = self._process_raw_srim_data(raw_srim_filepath)
            # Save the newly processed data
            if e_kev is not None and length_um is not None:
                np.savetxt(processed_srim_filepath, np.column_stack((e_kev, dee_dx, den_dx, length_um)), header="Energy(keV)    dEe/dx(keV/micro_m)  dEn/dx(keV/micro_m)  x(micro_m)", fmt="%.6e")

        if e_kev is None or length_um is None:
            print(f"Warning: Could not load or process SRIM data for Z={ion_z}. Skipping.")
            return None

        # Create interpolation function and cache it
        unique_length, unique_indices = np.unique(length_um, return_index=True)
        unique_e = e_kev[unique_indices]

        if len(unique_e) < 2:
            print(f"Error: Not enough unique data points for interpolation for Z={ion_z}. Skipping.")
            return None

        e_to_x_func = InterpolatedUnivariateSpline(unique_e, unique_length, k=1)
        self._srim_cache[ion_z] = (e_to_x_func, unique_e, dee_dx, den_dx, unique_length)
        return e_to_x_func, unique_e, dee_dx, den_dx, unique_length
    
    def _process_raw_srim_data(self, raw_srim_filepath):
        """
        Helper method to process raw SRIM data from a file.
        Returns (e_kev, length_nm) or (None, None) if file not found/error.
        """
        if not os.path.exists(raw_srim_filepath):
            print(f"Error: Raw SRIM file not found at {raw_srim_filepath}.")
            return None, None

        data = np.loadtxt(raw_srim_filepath, usecols=(0, 2, 3, 4, 6, 8), unpack=True)
        e_raw, dee_dx, den_dx, x_raw, y_raw, z_raw = data
        unit_e, unit_x = np.loadtxt(raw_srim_filepath, usecols=(1, 5), dtype=str, unpack=True)

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
                x_um[j], y_um[j], z_um[j] = x_raw[j], y_raw[j], z_raw[j]

        # Calculate total projected length
        length_um = np.sqrt(x_um**2 + y_um**2 + z_um**2)

        return e_kev, dee_dx, den_dx, length_um
    
    def calculate_fission_spectrum(self, x_bins):
        """
        Calculates the differential track rate from U-238 spontaneous fission.
        """
        print("Calculating spontaneous fission background...")
        
        Z_fission, A_fission, _ = self._load_nuclear_data("U238.dat")
        Z_bind, A_bind, B_bind = self._load_nuclear_data("BindingEne.txt")
        binding_map = {(int(z), int(a)): b for z, a, b in zip(Z_bind, A_bind, B_bind)}

        B0_U238 = 7.570126
        tau_U238_fission = 6.45e3  # Myr
        mass_U238 = 238.0
        u_concentration = self.config["uranium_concentration_g_g"]

        fission_rate_factor = (u_concentration * (AVOGADRO_NUMBER / mass_U238) * 1e3) / tau_U238_fission
        
        total_track_lengths_um = []
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
                track1 = srim_func1(Ek1_MeV * 1e3)
                track2 = srim_func2(Ek2_MeV * 1e3)
                total_track_lengths_um.append(track1 + track2)
                
        counts, bin_edges = np.histogram(total_track_lengths_um, bins=x_bins)
        bin_widths = np.diff(bin_edges)
        
        non_zero_widths = np.where(bin_widths > 0, bin_widths, 1)
        dRdx = (counts / num_events) * fission_rate_factor / non_zero_widths
        
        return dRdx

    def _calculate_flux_weights(self, scenario, energy_bins_gev, total_simulated_muons):
        """Integrates muon flux models over energy bins to get event weights."""
        
        scenario_name, flux_filename = scenario
        flux_filepath = os.path.join(self.data_path, "flux_data", f"{flux_filename}.txt")
        energy_GeV, flux_val = np.loadtxt(flux_filepath, usecols=(0, 1), unpack=True)                   
        flux_interpolator = log_interp1d(energy_GeV, flux_val / (energy_GeV**2))

        weights = []
        for i in range(len(energy_bins_gev) - 1):
            integrated_flux, _ = quad(flux_interpolator, energy_bins_gev[i], energy_bins_gev[i+1])
            weight = (integrated_flux * SECONDS_PER_MYR) / total_simulated_muons
            weights.append(weight)
        return weights

    def _get_all_fragments(self, energy_names_gev):
        """Dynamically determines the list of all nuclear fragments from Geant4 output."""
        geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_QGSP")
        nuclei_dir = os.path.join(geant4_input_dir, "Nuclei")
        all_fragments = set()
        
        for energy_name in energy_names_gev:
            filepath = os.path.join(nuclei_dir, f"outNuclei_{energy_name}.txt")
            if os.path.exists(filepath):
                names = np.loadtxt(filepath, usecols=0, dtype=str)
                for name in names:
                    clean_name = ''.join([char for char in name if char.isalpha()])
                    if clean_name not in ['He', 'alpha', 'proton']:
                        all_fragments.add(clean_name)
        return sorted(list(all_fragments))

    def _data_counter(self, er_bins, recoil_energies, weight):
        """Helper function to histogram weighted recoil energies."""
        counts, _ = np.histogram(recoil_energies, bins=er_bins)
        return counts * weight

    def _process_geant4_data(self, scenario, energy_names_gev, energy_bins_gev, target_thickness_cm, total_simulated_muons=1e4):
        """
        Processes raw Geant4 data for all scenarios and saves the results.
        This method replaces the main logic of the original 'Halite.py' script.
        """
        print("Processing raw Geant4 data for all scenarios...")
        
        weights = self._calculate_flux_weights(scenario, energy_bins_gev, total_simulated_muons)
        all_fragments = self._get_all_fragments(energy_names_gev)
        
        geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_QGSP")
        
        scenario_name, flux_filepath = scenario
        
        all_recoil_spectra = {}

        for nucleus_type in self.config["target_nuclei_geant4"]:
            nucleus_dir = os.path.join(geant4_input_dir, nucleus_type)

            if nucleus_type == "Nuclei":
                fragment_spectra = {frag: np.zeros(len(RECOIL_ENERGY_BINS_MEV) - 1) for frag in all_fragments}
                for i, energy_name in enumerate(energy_names_gev):
                    filepath = os.path.join(nucleus_dir, f"outNuclei_{energy_name}.txt")
                    if not os.path.exists(filepath): continue

                    names, _, energies = np.loadtxt(filepath, usecols=(0, 1, 2), dtype=str, unpack=True)
                    energies = energies.astype(float)

                    for name, energy in zip(names, energies):
                        clean_name = ''.join([char for char in name if char.isalpha()])
                        if clean_name in fragment_spectra:
                            bin_index = np.digitize(energy, RECOIL_ENERGY_BINS_MEV) - 1
                            if 0 <= bin_index < len(fragment_spectra[clean_name]):
                                fragment_spectra[clean_name][bin_index] += weights[i]
                all_recoil_spectra.update(fragment_spectra)
            else:
                total_spectrum = np.zeros(len(RECOIL_ENERGY_BINS_MEV) - 1)
                for i, energy_name in enumerate(energy_names_gev):
                    filepath = os.path.join(nucleus_dir, f"out{nucleus_type}_{energy_name}.txt")
                    if not os.path.exists(filepath): continue

                    recoil_energies_mev = np.loadtxt(filepath, usecols=2)
                    counts, _ = np.histogram(recoil_energies_mev, bins=RECOIL_ENERGY_BINS_MEV)
                    total_spectrum += counts * weights[i]
                all_recoil_spectra[nucleus_type] = total_spectrum

        # Normalize and save
        output_dir = os.path.join(self.data_path, "processed_recoils")
        os.makedirs(output_dir, exist_ok=True)
        output_filepath = os.path.join(output_dir, f"{self.name}_muon_recoil_{flux_filepath}.npz")

        bin_widths_mev = np.diff(RECOIL_ENERGY_BINS_MEV)
        norm_factor = (bin_widths_mev * target_thickness_cm * self.config['density_g_cm3'] * 1e-3)

        normalized_spectra = {}
        for name, spectrum in all_recoil_spectra.items():
            normalized_spectra[name] = np.divide(spectrum, norm_factor, out=np.zeros_like(spectrum), where=norm_factor!=0)

        np.savez(output_filepath, Er_bins=RECOIL_ENERGY_BINS_MEV, **normalized_spectra)
        print(f"    - Saved processed data to {output_filepath}")

    def _convert_recoil_to_track_spectrum(self, x_bins, recoil_data, depth_mwe):
        """Helper to convert a full dR/dEr spectrum to dR/dx."""
        attenuation = np.exp(-depth_mwe / 2500.0) # Attenuation in m.w.e.
        er_bins = recoil_data['Er_bins']
        er_mid_mev = er_bins[:-1] + np.diff(er_bins) / 2.0
        
        dRdx_by_nucleus = {}

        dRdx_total = np.zeros(len(x_bins) - 1)
        x_mid_nm = x_bins[:-1] + np.diff(x_bins) / 2.0
        
        # Loop through all fragment types defined for the mineral
        for nucleus_name in self.config['fragments']+self.config['target_nuclei_geant4']:
            if nucleus_name not in recoil_data: continue
            
            dRdEr_mev = recoil_data[nucleus_name]
            dRdEr_interp = interp1d(er_mid_mev, dRdEr_mev, bounds_error=False, fill_value=0.0)
            e_attenuated = attenuation*(er_mid_mev+510)-510
            dRdEr_atten = interp1d(er_mid_mev, dRdEr_interp(e_attenuated)*attenuation, bounds_error=False, fill_value=0.0)

            clean_name = ''.join([char for char in nucleus_name if char.isalpha()])

            ion_z = element(clean_name).atomic_number
            srim_func, e, dee_dx, den_dx, x = self._load_srim_data(ion_z)
            if not srim_func: continue

            sorted_indices = np.argsort(x)
            
            x_to_e_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, e[sorted_indices]*1e-3, k=1)
            x_to_dedx_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, (dee_dx[sorted_indices]+den_dx[sorted_indices])*1e-6, k=1)
            
            e_at_x = x_to_e_func(x_mid_nm)
            
            dRdx_nucleus = dRdEr_atten(e_at_x) * x_to_dedx_func(x_mid_nm)
            dRdx_by_nucleus[nucleus_name] = dRdx_nucleus
            dRdx_total += dRdx_nucleus
            
        dRdx_by_nucleus["total"] = dRdx_total
        
        return dRdx_by_nucleus

    def calculate_muon_signal_spectrum(self, x_bins, scenario, energy_names_gev, energy_bins_gev, target_thickness_cm, depth_mwe, total_simulated_muons=1e4, nucleus="total"):
        """
        Calculates the final muon-induced track length spectrum, including deposition model.
        """
        scenario_name, flux_filename = scenario
        
        if scenario_name in self._recoil_cache:
            recoil_data = self._recoil_cache[scenario_name]
        else:
            filepath = os.path.join(self.data_path, "processed_recoils", f"{self.name}_muon_recoil_{flux_filename}.npz")
            if not os.path.exists(filepath):
                self._process_geant4_data(scenario, energy_names_gev, energy_bins_gev, target_thickness_cm, total_simulated_muons)

            recoil_data = np.load(filepath)
            self._recoil_cache[scenario_name] = recoil_data

        dRdx_at_depth = self._convert_recoil_to_track_spectrum(x_bins, recoil_data, depth_mwe)
        
        if nucleus=="total":
            return dRdx_at_depth["total"]
        elif nucleus=="all":
            return dRdx_at_depth
        else:
            return dRdx_at_depth[nucleus]

    
    def integrate_muon_signal_spectrum(self, x_bins, scenario, energy_names_gev, energy_bins_gev, exposure_window_myr, sample_mass_kg, target_thickness_cm, initial_depth = 0, deposition_rate_m_kyr=0, overburden_density_g_cm3=1., total_simulated_muons=1e4):
        
        time_steps_myr = np.linspace(0, exposure_window_myr, 20)
        dt_myr = time_steps_myr[1] - time_steps_myr[0]
        total_tracks = np.zeros(len(x_bins) - 1)
     
        for t_myr in time_steps_myr:
            print(f"  - Processing timestep: {t_myr*1000} kyr")

            depth_mwe = initial_depth + deposition_rate_m_kyr * (t_myr * 1e3) * overburden_density_g_cm3
            
            dRdx_at_depth = self.calculate_muon_signal_spectrum(x_bins, scenario, energy_names_gev, energy_bins_gev, target_thickness_cm, depth_mwe, total_simulated_muons)
            
            total_tracks += np.random.poisson(dRdx_at_depth * sample_mass_kg * dt_myr * np.diff(x_bins))

        return total_tracks