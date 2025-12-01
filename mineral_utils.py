import numpy as np
import os
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.integrate import quad
from mendeleev import element
from tqdm import tqdm
from multiprocessing import Pool

# --- Physical Constants ---
PROTON_MASS_MEV = 938.3
NEUTRON_MASS_MEV = 939.6
AVOGADRO_NUMBER = 6.022e23
MYR_PER_SECOND = 1/ (60 * 60 * 24 * 365 * 1e6)

# --- Binning setup ---

RECOIL_N_BINS= 101
RECOIL_ER_MIN_LOG_MEV= -2 
RECOIL_ER_MAX_LOG_MEV= 3
RECOIL_ENERGY_BINS_MEV = np.logspace(RECOIL_ER_MIN_LOG_MEV, RECOIL_ER_MAX_LOG_MEV, RECOIL_N_BINS)

LENGTH_N_BINS = 1000
LENGTH_MIN_LOG_NM = 1.5
LENGTH_MAX_LOG_NM = 5.5
TRACK_LENGTH_BINS_NM = np.logspace(LENGTH_MIN_LOG_NM, LENGTH_MAX_LOG_NM, LENGTH_N_BINS)

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

def _asymmetric_gaussian_kernel(size, sigma_left, sigma_right=None):
    """
    Generates an asymmetric Gaussian kernel for convolution.

    If sigma_left and sigma_right are not provided, it generates a standard
    symmetric Gaussian kernel.

    Args:
        size (int): The total size (number of points) of the kernel. Must be odd.
        sigma_left (float): The standard deviation for the left tail.
        sigma_right (float, optional): The standard deviation for the right tail. If None, the kernel is a symmetric gaussian with sigma = sigma_left.

    Returns:
        np.ndarray: The normalized 1D convolution kernel.
    """
    if size % 2 == 0:
        raise ValueError("Kernel size must be odd.")

    center = size // 2
    x = np.arange(size)

    if sigma_right is None:
        sigma_right = sigma_left
    
    kernel = np.zeros(size)
    
    kernel[:center] = np.exp(-(x[:center] - center)**2 / (2 * sigma_left**2))
    kernel[center] = 1.0
    kernel[center+1:] = np.exp(-(x[center+1:] - center)**2 / (2 * sigma_right**2))
    
    return kernel / np.sum(kernel)

# --- Main Paleodetector Class ---

class Paleodetector:
    """
    A class to handle all physics calculations for a specific mineral.
    
    This class is initialized with a configuration dictionary and provides
    methods to calculate various signal and background components for
    paleo-detector analysis.
    """
    
    def __init__(self, mineral_config, data_path_prefix="Data"):
        """
        Initializes the Paleodetector object and sets up its properties and data caches.

        Args:
            mineral_config (dict): A dictionary containing all properties of the mineral,
                                   such as name, composition, and nuclear data.
            data_path_prefix (str, optional): The relative path to the main data directory.
                                              Defaults to "Data/".
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
        self._neutron_bkg_cache = {}
        self._flux_interpolators = {}
        self._energy_GeV = {}
        self._depth_interpolators = {}
        
        print(f"Initialized Paleodetector: {self.name}")

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
    
    def _load_neutron_bkg(self):
        """
        Loads and caches the radiogenic neutron background spectrum from a file.

        Returns:
            dict: A dictionary of interpolation functions, keyed by nucleus symbol,
                  for the differential neutron recoil rate [events/kg/day/keV].
        """
        if self._neutron_bkg_cache:
            return self._neutron_bkg_cache

        fname = os.path.join(self.data_path, "neutron_data", f"{self.name}_ninduced_wan.dat")
        if not os.path.exists(fname):
            raise FileNotFoundError(f"Neutron background file not found: {fname}")

        with open(fname) as f:
            head = f.readlines()[1]
            columns = [c.strip() for c in head.split(",")]
        
        data = np.loadtxt(fname)
        E_list = data[:, 0]
        
        neutron_interp_dict = {}
        for i, nuc in enumerate(self.config['nuclei']):
            dRdE_list = np.zeros_like(E_list)
            for j, col_name in enumerate(columns):
                if col_name.startswith(nuc):
                    dRdE_list += data[:, j]
            neutron_interp_dict[nuc] = interp1d(E_list, dRdE_list, bounds_error=False, fill_value=0.0)
        
        self._neutron_bkg_cache = neutron_interp_dict
        return self._neutron_bkg_cache
    
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
            print(f"Pre-processed SRIM data for {ion_symbol} (Z={ion_z}) not found. Processing raw file.")
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

        e_to_x_func = InterpolatedUnivariateSpline(unique_e, unique_length, k=1)
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

        data = np.genfromtxt(raw_srim_filepath, usecols=(0, 2, 3, 4, 6, 8), unpack=True, skip_header=26, skip_footer=13, converters={0: comma_to_dot, 2: comma_to_dot, 3: comma_to_dot, 4: comma_to_dot, 6: comma_to_dot, 8: comma_to_dot})

        e_raw, dee_dx, den_dx, x_raw, y_raw, z_raw = data
        unit_e, unit_x = np.genfromtxt(raw_srim_filepath, dtype=str, skip_header=26, skip_footer=13, usecols=(1, 5), unpack=True)

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

    def calculate_neutron_spectrum(self, x_bins=TRACK_LENGTH_BINS_NM):
        """
        Calculates the differential track rate from radiogenic neutrons.

        Args:
            x_bins (np.ndarray): The bin edges for the output track length spectrum [nm].

        Returns:
            np.ndarray: The differential track rate (dR/dx) [events/kg/Myr/nm].
        """
        if not self.config.get("uranium_concentration_g_g"):
            print("Warning: Uranium concentration not set. Cannot calculate neutron background.")
            return np.zeros(len(x_bins) - 1)
            
        print("Calculating radiogenic neutron background...")
        neutron_interp_dict = self._load_neutron_bkg()
        x_mid = x_bins[:-1] + np.diff(x_bins) / 2.0
        dRdx = np.zeros_like(x_mid)

        for i, nuc_name in enumerate(self.config['nuclei']):
            if nuc_name == "H":
                continue

            ion_z = element(nuc_name).atomic_number
            srim_func, e, dee_dx, den_dx, x = self._load_srim_data(ion_z)
            if not srim_func: continue

            sorted_indices = np.argsort(x)
        
            x_to_e_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, e[sorted_indices], k=1)
            x_to_dedx_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, (dee_dx[sorted_indices]+den_dx[sorted_indices])*1e-3, k=1)

            dRdE_kev = neutron_interp_dict[nuc_name](x_to_e_func(x_mid))
            dRdx += self.config['stoich'][i] *dRdE_kev * np.abs(x_to_dedx_func(x_mid))
    
        return dRdx * self.config["uranium_concentration_g_g"] / 0.1e-9
    
    def integrate_neutron_spectrum(self, x_bins, age, sample_mass, x_grid=TRACK_LENGTH_BINS_NM):

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
        x_mids_grid = x_grid[:-1] + np.diff(x_grid) / 2.0

        drdx = self.calculate_neutron_spectrum(x_grid) * age * sample_mass

        total_tracks_interp  = interp1d(x_mids_grid, np.array(drdx),  bounds_error=False, fill_value='extrapolate')

        total_tracks = [quad(total_tracks_interp, x_bins[i], x_bins[i+1])[0] for i in range(len(x_mids))]

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
            
                x_to_e_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, e[sorted_indices], k=1)
                x_to_dedx_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, (dee_dx[sorted_indices]+den_dx[sorted_indices])*1e-3, k=1)

                dRdE_kev = np.vectorize(DMU.dRdE_CEvNS)(x_to_e_func(x_mid), element(nuc_name).protons, element(nuc_name).neutrons, flux_name=flux_name)

                dRdx += self.config['stoich'][i] * dRdE_kev * np.abs(x_to_dedx_func(x_mid))
        
        return dRdx * 365 * 1e6
    
    def integrate_nu_spectrum(self, x_bins, age, sample_mass, flux_name="all", x_grid=TRACK_LENGTH_BINS_NM):

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
        x_mids_grid = x_grid[:-1] + np.diff(x_grid) / 2.0

        drdx = self.calculate_nu_spectrum(x_grid, flux_name) * age * sample_mass

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
    
    def integrate_fission_spectrum(self, x_bins, age, sample_mass, x_grid=TRACK_LENGTH_BINS_NM):

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
        x_mids_grid = x_grid[:-1] + np.diff(x_grid) / 2.0

        drdx = self.calculate_fission_spectrum(x_grid) * age * sample_mass

        total_tracks_interp  = interp1d(x_mids_grid, np.array(drdx),  bounds_error=False, fill_value='extrapolate')

        total_tracks = [quad(total_tracks_interp, x_bins[i], x_bins[i+1])[0] for i in range(len(x_mids))]

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
        pri_energies = np.logspace(-3, 3, 10)

        tab_species = 'mu-' if species == 'mu-' or species == 'mu+' else 'neutron'

        for energy_name in pri_energies:

            filepath = os.path.join(self.data_path, "Geant4_data", f"StdRock_{tab_species}", f"outNuclei_{energy_name:.6f}.txt")
            if not os.path.exists(filepath):
                print(f"Geant4 data file not found: {filepath}")
                continue

            energies, depth, rem_energy = np.loadtxt(filepath, usecols=(2, 3, 5), dtype = str, unpack=True)
            energies = energies.astype(float)
            depth = 500 - depth.astype(float)/1000.
            rem_energy = rem_energy.astype(float)

            pen.append(depth[rem_energy == 0].mean()*2.65)
            width.append(depth[rem_energy == 0].std()*2.65)

        self._depth_interpolators[species] = {}
        maxdepth = interp1d(pri_energies[:-1], pen[:-1], kind='linear', fill_value='extrapolate', bounds_error=False)
        meanwidth = interp1d(pri_energies[:-1], width[:-1], kind='linear', fill_value='extrapolate', bounds_error=False)
    
        self._depth_interpolators[species]['maxdepth'] = lambda x : np.clip(maxdepth(x), 0.5e-3, np.inf)
        self._depth_interpolators[species]['meanwidth']= lambda x : np.clip(meanwidth(x), 1.e-4, np.inf)

    def _interpolate_flux_scenarios(self, scenario_config, species='mu-'):
        """
        Interpolates particle flux data for a given scenario configuration.
        
        Args:
            scenario_config (dict): A dictionary containing the scenario configuration,
                                    including 'name' and 'event_fluxes'.
            species (str, optional): The particle species to simulate ('mu+', 'mu-', or 'neutron'). Defaults to 'mu-'.

        """
        times = []
        flux_arrays = []

        if species == 'mu+' or species == 'mu-':
            col = 1
        elif species == 'neutron':
            col = 2
        else:
            raise ValueError('Accepted species are mu+, mu- and neutron')
        
        for scenario in scenario_config['event_fluxes'].items():
            times.append(scenario[0])
            _, filename_tag = scenario[1]
            flux_filepath = os.path.join(self.data_path, "flux_data", f"{filename_tag}.txt")
            if not os.path.exists(flux_filepath):
                print(f"Flux file not found: {flux_filepath}")
                continue
            
            energies, flux = np.loadtxt(flux_filepath, usecols=(0, col), unpack=True)
            if species == 'mu+' or species == 'mu-':
                flux /= 2.
            flux_arrays.append(flux)
        
        flux_arrays = np.array(flux_arrays)
        if len(times)== 1:
            kind = "nearest"
        else:
            kind = "linear"

        interpolators = []
        for i in range(len(flux_arrays[0])):
            interp_func = interp1d(times, flux_arrays[:, i], kind=kind, fill_value="extrapolate")
            interpolators.append(interp_func)
        
        self._flux_interpolators[f'{scenario_config["name"]}_{species}'] = interpolators
        self._energy_GeV[scenario_config["name"]] = energies

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
            if os.path.exists(filepath):
                names = np.loadtxt(filepath, usecols=0, dtype=str)
                for name in names:
                    if name not in ['He3', 'He4', 'He5', 'He6', 'He7', 'He8', 'alpha', 'proton', 'neutron']:
                        all_fragments.add(name)
        return sorted(list(all_fragments))

    def _process_geant4_data(self, t_kyr, scenario_name, energy_bins_gev, depth_mwe=0., total_simulated_particles=1e4, target_thickness_mm=5., species='mu-'):
        """
        Processes raw Geant4 data for a given scenario, creating a normalized recoil spectrum file.

        Args:
            t_kyr (float): The time in kiloyears for which to process the data.
            scenario_name (str): The name of the flux scenario to use.
            energy_bins_gev (np.ndarray): The energy bin edges [GeV].
            target_thickness_mm (float): The thickness of the target in the Geant4 simulation [mm].
            depth_mwe (float, optional): Shielding depth [m.w.e.]. Defaults to 0.
            total_simulated_particles (float, optional): Number of particles per Geant4 run. Defaults to 1e4.
            species (str, optional): The particle species to simulate ('mu+', 'mu-', or 'neutron'). Defaults to 'mu-'.
        """
        
        depth_min = depth_mwe
        depth_max = depth_mwe + target_thickness_mm * 0.001 * self.config['density_g_cm3']

        if not self._flux_interpolators[f'{scenario_name}_{species}']:
            raise ValueError(f"Flux interpolators not initialized for scenario {scenario_name} and species {species}.")

        if not self._depth_interpolators[species]:
            self._load_depth_interpolators(species)

        flux_val = np.array([interp_func(t_kyr) for interp_func in self._flux_interpolators[f'{scenario_name}_{species}']])

        flux_interpolator = log_interp1d(self._energy_GeV[scenario_name], flux_val)
        maxdepth = self._depth_interpolators[species]['maxdepth']
        meanwidth = self._depth_interpolators[species]['meanwidth']

        all_fragments = self._get_all_fragments(energy_bins_gev[:-1], species)

        geant4_input_dir = os.path.join(self.data_path, "Geant4_data", f"{self.name}_{species}")
                
        all_recoil_spectra = {}

        fragment_spectra = {frag: np.zeros(len(RECOIL_ENERGY_BINS_MEV) - 1) for frag in all_fragments}
    
        for i in range(len(energy_bins_gev) - 1):
            e_min = energy_bins_gev[i]
            e_max = energy_bins_gev[i+1]
            weight_elastic = quad(lambda e: quad(lambda x: np.clip(np.exp(-(x - maxdepth(e))/meanwidth(e)), 0. , 1.) / (maxdepth(e)+meanwidth(e)), depth_min, depth_max)[0]*flux_interpolator(e), e_min, e_max)[0]
            
            if species == 'mu-':
                weight_peak = quad(lambda e: quad(lambda x: np.exp(-(x-maxdepth(e))**2/(2*meanwidth(e)**2))/(np.sqrt(2*np.pi)*meanwidth(e)), depth_min, depth_max)[0]*flux_interpolator(e), e_min, e_max)[0]

            filepath = os.path.join(geant4_input_dir, f"outNuclei_{e_min:.6f}.txt")
            if not os.path.exists(filepath): continue

            names, rec_energies, rem_energies = np.loadtxt(filepath, usecols=(0, 2, 5), dtype=str, unpack=True)
            rec_energies = rec_energies.astype(float)
            rem_energies = rem_energies.astype(float)

            for name, rec_energy, rem_energy in zip(names, rec_energies, rem_energies):
                if name in fragment_spectra:
                    bin_index = np.digitize(rec_energy, RECOIL_ENERGY_BINS_MEV) - 1
                    if 0 <= bin_index < len(fragment_spectra[name]):
                        if rem_energy == 0.0 and species == 'mu-':
                            fragment_spectra[name][bin_index] += weight_peak
                        else:
                            fragment_spectra[name][bin_index] += weight_elastic
        all_recoil_spectra.update(fragment_spectra)

        output_dir = os.path.join(self.data_path, "processed_recoils")
        os.makedirs(output_dir, exist_ok=True)
        output_filepath = os.path.join(output_dir, f"{self.name}_{species}_recoil_{scenario_name}_{t_kyr}kyr_{depth_mwe:.1f}mwe.npz")

        bin_widths_mev = np.diff(RECOIL_ENERGY_BINS_MEV)
        norm_factor = (bin_widths_mev * target_thickness_mm * self.config['density_g_cm3'] * total_simulated_particles * MYR_PER_SECOND)

        normalized_spectra = {}
        for name, spectrum in all_recoil_spectra.items():
            normalized_spectra[name] = np.divide(spectrum, norm_factor, out=np.zeros_like(spectrum), where=norm_factor!=0)

        np.savez(output_filepath, Er_bins=RECOIL_ENERGY_BINS_MEV, **normalized_spectra)
        print(f"    - Saved processed data to {output_filepath}")

    def _convert_recoil_to_track_spectrum(self, x_bins, recoil_data, energy_bins_gev, species='mu-'):
        """
        Converts a full differential recoil energy spectrum (dR/dEr) to a track length spectrum (dR/dx).

        Args:
            x_bins (np.ndarray): The bin edges for the output track length spectrum [nm].
            recoil_data (np.lib.npyio.NpzFile): The loaded .npz file containing recoil energy spectra.
            energy_bins_gev (np.ndarray): The energy bin edges [GeV].
            species (str, optional): The particle species to simulate ('mu+', 'mu-', or 'neutron'). Defaults to 'mu-'.

            
        Returns:
            dict: A dictionary of differential track rates (dR/dx) [events/kg/Myr/nm],
                  keyed by nucleus/fragment symbol, including a "total" key.
        """
        er_bins = recoil_data['Er_bins']
        er_mid_mev = er_bins[:-1] + np.diff(er_bins) / 2.0
        all_fragments = self._get_all_fragments(energy_bins_gev[:-1], species)

        dRdx_by_nucleus = {}

        dRdx_total = np.zeros(len(x_bins) - 1)
        x_mid_nm = x_bins[:-1] + np.diff(x_bins) / 2.0
        
        for nuclide_name in all_fragments:
            if nuclide_name not in recoil_data: continue
            
            dRdEr_mev = recoil_data[nuclide_name]
            dRdEr_interp = interp1d(er_mid_mev, dRdEr_mev, bounds_error=False, fill_value=0.0)

            nucleus_name = ''.join(filter(str.isalpha, nuclide_name))

            ion_z = element(nucleus_name).atomic_number
            srim_func, e, dee_dx, den_dx, x = self._load_srim_data(ion_z)
            if not srim_func: continue

            sorted_indices = np.argsort(x)
            
            x_to_e_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, e[sorted_indices]*1e-3, k=1)
            x_to_dedx_func = InterpolatedUnivariateSpline(x[sorted_indices]*1e3, (dee_dx[sorted_indices]+den_dx[sorted_indices])*1e-6, k=1)
            
            e_at_x = x_to_e_func(x_mid_nm)
            
            dRdx_nucleus = dRdEr_interp(e_at_x) * x_to_dedx_func(x_mid_nm)
            dRdx_by_nucleus[nuclide_name] = dRdx_nucleus
            dRdx_total += dRdx_nucleus
            
        dRdx_by_nucleus["total"] = dRdx_total
        
        return dRdx_by_nucleus

    def calculate_particle_signal_spectrum(self, x_bins, t_kyr, scenario_name, energy_bins_gev, depth_mwe, total_simulated_particles=1e5,  target_thickness_mm=5., species='mu-', nucleus="total", time_precision=0):
        """
        Calculates the final particle-induced differential track length spectrum (dR/dx) for a given depth.

        Args:
            x_bins (np.ndarray): The bin edges for the output track length spectrum [nm].
            t_kyr (float): The time in kiloyears for which to calculate the spectrum.
            scenario_name (str): The name of the flux scenario to use.
            energy_bins_gev (np.ndarray): The energy bin edges [GeV].
            target_thickness_mm (float): Target thickness in the Geant4 simulation [mm].
            depth_mwe (float): Shielding depth [m.w.e.].
            total_simulated_particles (float, optional): Number of particles per Geant4 run. Defaults to 1e5.
            nucleus (str, optional): Which nucleus/fragment spectrum to return ('total', 'all', or a specific symbol). Defaults to "total".
            time_precision (int, optional): Decimal precision for time-based filenames. Defaults to 0 (O(kyr)).
            species (str, optional): The particle species to simulate ('mu+', 'mu-', or 'neutron'). Defaults to 'mu-'.

            
        Returns:
            np.ndarray or dict: The differential track rate(s) (dR/dx) [events/kg/Myr/nm].
        """
        t_kyr = round(t_kyr, time_precision)

        filepath = os.path.join(self.data_path, "processed_recoils", f"{self.name}_{species}_recoil_{scenario_name}_{t_kyr}kyr_{depth_mwe:.1f}mwe.npz")
        if not os.path.exists(filepath):
            self._process_geant4_data(t_kyr, scenario_name, energy_bins_gev, depth_mwe, total_simulated_particles, target_thickness_mm, species)

        recoil_data = np.load(filepath)
        self._recoil_cache[f'{scenario_name}_{species}'] = recoil_data

        dRdx_at_depth = self._convert_recoil_to_track_spectrum(x_bins, recoil_data, energy_bins_gev, species)
        
        if nucleus=="total":
            return dRdx_at_depth["total"]
        elif nucleus=="all":
            return dRdx_at_depth
        else:
            return dRdx_at_depth[nucleus]
    
    def _integration_worker(self, args):
        """
        Helper worker function for parallel processing that calculates track counts for a single timestep.

        Args:
            args (tuple): A tuple containing all necessary arguments for a single timestep calculation.

        Returns:
            np.ndarray: The number of tracks produced in each bin for this timestep.
        """
        x_bins, t_kyr, scenario_name, energy_bins_gev, \
        initial_depth, deposition_rate_m_kyr, \
        overburden_density_g_cm3, total_simulated_particles, target_thickness_mm, species = args

        depth_mwe = initial_depth + deposition_rate_m_kyr * t_kyr * overburden_density_g_cm3

        dRdx_at_depth = self.calculate_particle_signal_spectrum(
            x_bins, t_kyr, scenario_name, energy_bins_gev,
            depth_mwe, total_simulated_particles, target_thickness_mm, species
        )

        return dRdx_at_depth, t_kyr

    def integrate_particle_signal_spectrum_parallel(
            self, 
            x_bins, 
            scenario_config, 
            energy_bins_gev, 
            exposure_window_kyr, 
            sample_mass_kg, 
            initial_depth=0, 
            deposition_rate_m_kyr=0, 
            overburden_density_g_cm3=1., 
            nsteps=None, 
            total_simulated_particles=1e5, 
            target_thickness_mm=5., 
            x_grid=TRACK_LENGTH_BINS_NM, 
            species='mu-'):
        """
        Calculates the final particle-induced track length spectrum by parallelizing the time integration.

        Args:
            x_bins (np.ndarray): The bin edges for the output track length spectrum [nm].
            scenario_config (dict): Configuration dictionary for the flux scenario.
            energy_bins_gev (np.ndarray): The energy bin edges [GeV].
            exposure_window_kyr (float): The total exposure time in kiloyears.
            sample_mass_kg (float): The mass of the sample in kilograms.
            initial_depth (float, optional): Initial depth in meters water equivalent [m.w.e.]. Defaults to 0.
            deposition_rate_m_kyr (float, optional): Deposition rate in meters per kiloyear. Defaults to 0.
            overburden_density_g_cm3 (float, optional): Overburden density in g/cm³. Defaults to 1.
            nsteps (int, optional): Number of time steps for integration. Defaults to 75*(number of flux changes in scenario_config).
            total_simulated_particles (float, optional): Number of particles per Geant4 run. Defaults to 1e4.        
            target_thickness_mm (float, optional): Thickness of the target [mm]. Defaults to 5.
            x_grid (np.ndarray, optional): The bin edges for the internal track length spectrum [nm]. Defaults to TRACK_LENGTH_BINS_NM.
            species (str, optional): The particle species to simulate ('mu+', 'mu-', or 'neutron'). Defaults to 'mu-'.
            
        Returns:
            np.ndarray: The total number of tracks expected in each track length bin.
        """

        self._interpolate_flux_scenarios(scenario_config, species)
        self._load_depth_interpolators(species)

        if not nsteps:
            nsteps = 75 * len(scenario_config["event_fluxes"])

        time_bins_kyr = np.linspace(0., exposure_window_kyr, nsteps + 1)

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
        x_mids_grid = x_grid[:-1] + np.diff(x_grid) / 2.0

        tasks = [(x_grid, t_kyr, scenario_config["name"], energy_bins_gev, 
                   initial_depth, deposition_rate_m_kyr, 
                  overburden_density_g_cm3, total_simulated_particles, target_thickness_mm, species)
                 for t_kyr in time_bins_kyr]

        with Pool() as pool:
            results = list(tqdm(pool.imap(self._integration_worker, tasks), total=len(tasks)))
        
        drdx_array, t_kyr = zip(*results)

        drdx_array = np.array(drdx_array)
        if len(time_bins_kyr) == 1:
            kind = "nearest"
        else:
            kind = "linear"

        drdx_interpolators = [interp1d(t_kyr, drdx, bounds_error=False, fill_value=0.0, kind=kind) for drdx in drdx_array.T]

        total_drdx = [quad(drdx_interpolator, 0, exposure_window_kyr)[0] for drdx_interpolator in drdx_interpolators]

        total_tracks_interp  = interp1d(x_mids_grid, np.array(total_drdx) * sample_mass_kg * 1e-3 * np.pi,  bounds_error=False, fill_value='extrapolate')

        total_tracks = [quad(total_tracks_interp, x_bins[i], x_bins[i+1])[0] for i in range(len(x_mids))]

        return total_tracks
    
def smear_spectrum(counts, size, sigma_left, sigma_right=None):
    """
    Applies asymmetric gaussian smearing to the track length distribution.

    Args:
        counts (np.ndarray): Track counts in the bins.
        size (int): The total size (number of points) of the kernel. Must be odd.
        sigma_left (float): The standard deviation for the left tail.
        sigma_right (float, optional): The standard deviation for the right tail. If None, the kernel is a symmetric gaussian with sigma = sigma_left.

    Returns:
        np.ndarray: The smeared track counts.
    """        
    smeared_counts = np.convolve(counts, _asymmetric_gaussian_kernel(size, sigma_left, sigma_right), mode='same')

    return smeared_counts

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

def slice_spectrum(x_bins, counts, angular_pdf=None, phi_cut_deg=0., l_min_measurable=300., l_max_measurable=50000., pit_width=500., bulk_etching_depth=100., f_phi= lambda phi: 1., n_samples=1e6, correction=True):
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
        l_min_measurable (float): Minimum measurable segment length, L_min [nm]. Tracks shorter than 
                                    this are lost due resolution limits.
        l_max_measurable (float): Maximum measurable segment length, L_max [nm]. This caps 
                                    the pit size due to saturation effects in the etching process.
        pit_width (float): Typical width of the etched pit [nm]. 
                            Tracks with parallel footprint smaller than this threshold will be measured by this.
        bulk_etching_depth (float): Minimum vertical development of the track [nm]. 
                                    Tracks shallower than this are lost due to etching away.
        f_phi (callable): Function applied to the segment length L_seg * f_phi(phi). Corrects 
                            for anisotropic enlargement (e.g., set to lambda phi: 1.0 for plasma etching).
        n_samples (int): Number of Monte Carlo tracks to simulate for accurate statistics.
        correction (bool): If True, applies pit_width correction. If not, the pit_width correction is ignored.

    Returns:
        np.ndarray: The resulting measured track count histogram N(L_meas), normalized to 
                    the total input counts.
    """
    x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0
    phi_cut_rad = np.deg2rad(phi_cut_deg)

    samples = np.random.choice(x_mids, size=int(n_samples), p=counts/np.sum(counts))
    
    phi_grid = np.linspace(0, np.pi / 2, 1000)

    if not angular_pdf:
        angular_pdf = np.sin(phi_grid)

    sampled_angles = np.random.choice(phi_grid, size=int(n_samples), p=angular_pdf / np.sum(angular_pdf))

    is_retained = sampled_angles >= phi_cut_rad

    samples_retained = samples[is_retained]
    phi_retained = sampled_angles[is_retained]

    seg_samples = np.random.uniform(low=0., high=samples_retained)

    measured_samples = seg_samples * f_phi(phi_retained)

    is_retained_length = measured_samples >= l_min_measurable

    measurable_samples = measured_samples[is_retained_length]
    measurable_angles = phi_retained[is_retained_length]

    is_retained_depth = measurable_samples * np.sin(measurable_angles) >= bulk_etching_depth

    if correction:
        corrected_measurable_samples = np.where(measurable_samples[is_retained_depth] * np.cos(measurable_angles[is_retained_depth]) >= pit_width/2., (pit_width/2.)+measurable_samples[is_retained_depth] * np.cos(measurable_angles[is_retained_depth]), pit_width)
    else:
        corrected_measurable_samples = measurable_samples[is_retained_depth]

    final_measured_samples = np.minimum(corrected_measurable_samples, l_max_measurable)

    hist_measured, _ = np.histogram(final_measured_samples, bins=x_bins, density=False)

    hist_norm = hist_measured * (np.sum(counts) / n_samples)

    return hist_norm