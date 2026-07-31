import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d

DEFAULT_ANGLES_DEG = np.rad2deg(np.arccos(np.linspace(1.0, 0.0, 11)))


class FluxTemplate:
    """
    A single precomputed (time-since-onset x energy) flux grid for one or
    more particle species, produced by one MCEq configuration.

    A template is either:
      * "steady" (`is_steady=True`): a single energy spectrum, constant in
        time. Used for baselines.
      * "transient": a 2D (t_since_kyr, energy_gev) grid representing the
        flux *excess* relative to a baseline, decaying back towards zero
        as t_since_kyr grows. Used for discrete events.

    Templates loaded from disk involve no MCEq / crflux import and are
    essentially instantaneous; only `get_or_build` (on a cache miss) pays
    the MCEq cost.
    """

    def __init__(self, name, energy_gev, species_grids, t_since_kyr=None):
        self.name = name
        self.energy_gev = np.asarray(energy_gev, dtype=float)
        self.t_since_kyr = None if t_since_kyr is None else np.asarray(t_since_kyr, dtype=float)
        self.is_steady = self.t_since_kyr is None

        self._grids = {}
        self._interp = {}
        for species, grid in species_grids.items():
            grid = np.asarray(grid, dtype=float)
            self._grids[species] = grid
            self._interp[species] = self._build_interpolator(grid)

    def _build_interpolator(self, grid):
        log_e = np.log10(self.energy_gev)
        log_flux = np.log10(np.clip(grid, 1e-300, None))
        if self.is_steady:
            return interp1d(
                log_e, log_flux, kind='linear',
                bounds_error=False, fill_value='extrapolate',
            )
        return RegularGridInterpolator(
            (self.t_since_kyr, log_e), log_flux,
            method='linear', bounds_error=False, fill_value=0.0,
        )

    @property
    def species(self):
        return list(self._grids.keys())

    def flux(self, species, energy_gev, t_since_kyr=None):
        """
        Evaluate this template for one species.

        Args:
            species: one of the species stored in this template.
            energy_gev: energies [GeV], array-like.
            t_since_kyr: required for transient templates (ignored for
                steady ones); broadcastable against energy_gev.
        """
        energy_gev = np.asarray(energy_gev, dtype=float)
        if species not in self._interp:
            return np.zeros_like(energy_gev)

        log_e = np.log10(energy_gev)
        if self.is_steady:
            return np.power(10.0, self._interp[species](log_e))

        if t_since_kyr is None:
            raise ValueError(f"Template '{self.name}' is transient: t_since_kyr is required.")

        t = np.asarray(t_since_kyr, dtype=float)

        t_clipped = np.clip(t, self.t_since_kyr[0], self.t_since_kyr[-1])
        t_b, e_b = np.broadcast_arrays(t_clipped, log_e)
        pts = np.stack([t_b.ravel(), e_b.ravel()], axis=-1)

        out = np.power(10.0, self._interp[species](pts)).reshape(t_b.shape)
        t_orig_b = np.broadcast_to(t, t_b.shape)
        return np.where(t_orig_b < 0.0, 0.0, out)

    def save(self, filepath):
        directory = os.path.dirname(filepath)
        if directory:
            os.makedirs(directory, exist_ok=True)
        payload = dict(
            energy_gev=self.energy_gev,
            is_steady=np.array(self.is_steady),
        )
        if self.t_since_kyr is not None:
            payload['t_since_kyr'] = self.t_since_kyr
        for species, grid in self._grids.items():
            payload[f'species__{species}'] = grid
        np.savez_compressed(filepath, **payload)

    @classmethod
    def load(cls, filepath):
        """Load a template previously written with `.save()`. No MCEq needed."""
        with np.load(filepath, allow_pickle=False) as data:
            energy_gev = data['energy_gev']
            t_since_kyr = data['t_since_kyr'] if 't_since_kyr' in data else None
            species_grids = {
                key.split('__', 1)[1]: data[key]
                for key in data.files if key.startswith('species__')
            }
        name = os.path.splitext(os.path.basename(filepath))[0]
        return cls(name, energy_gev, species_grids, t_since_kyr=t_since_kyr)

    @classmethod
    def get_or_build(cls, params, template_dir="Data/flux_data", force_rebuild=False):
        """
        `params` is a small dict describing a
        physical scenario, e.g. `{"kind": "SN", "distance_pc": 200}`.

        If it's already on disk (and `force_rebuild` is False), loads
        and returns it. Otherwise, calls MCEq (lazily imported) to build it, saves it
        to `template_dir`, and returns it. Next time, this scenario is loaded.
        """
        name = "_".join(str(val) for val in params.values())
        path = os.path.join(template_dir, name + ".npz")
        if os.path.exists(path) and not force_rebuild:
            return cls.load(path)
        params_tuple = tuple(params.values())[1:]
        if params.get("kind") == "Baseline":
            model = None
        elif params.get("kind") == "SN":
            model = SNHG12
        elif params.get("kind") == "Egal":
            model = EgalHG12
        elif params.get("kind") == "Enhanced":
            model = EnhancedHG12
        else:
            raise ValueError(f"Unknown template kind: {params.get('kind')}")

        tcut = params.get("tcut", 500.)
        t_since_kyr = np.linspace(0.001, tcut, 25) * 1.e3

        template = _build_template(name=name, model=model, params=params_tuple, t_since_kyr=t_since_kyr)
        template.save(path)
        return template

class FluxHistory:
    """
    Baseline template + a list of discrete events (plain dicts), summed on
    demand. Structurally the flux-side analogue of `overburden_history`:
    a small spec goes in, one object answers every downstream query.

    Timeline: see module docstring -- 0 = present, negative = past.
    """

    def __init__(self, baseline="Baseline", events=None, template_dir="Data/flux_data"):
        """
        Args:
            baseline: a FluxTemplate, a cached template name (str), or a
                params dict (built on demand) -- e.g. {"kind": "baseline"}.
            events: list of dicts, each with 'start_time_kyr' (or
                'time_kyr') and 'template' (FluxTemplate / name / params
                dict), and optionally 'weight' / 'label'.
            template_dir: directory used to load/save templates.
            name: optional human label; if omitted, `.signature` derives
                one from the actual content.
        """
        self.template_dir = template_dir
        self.baseline = self._resolve(baseline)
        self.events = []
        for event in (events or []):
            start_time = event.get("start_time_kyr")
            self.add_event(start_time, event["template"], weight=event.get("weight", 1.0), label=event.get("label"))

    def _resolve(self, template):
        if isinstance(template, FluxTemplate):
            return template
        if isinstance(template, dict):
            return FluxTemplate.get_or_build(template, self.template_dir)
        if isinstance(template, str):
            path = os.path.join(self.template_dir, f"{template}.npz")
            if os.path.exists(path):
                return FluxTemplate.load(path)
            raise FileNotFoundError(
                f"No cached template named '{template}' in {self.template_dir}. "
                f"Pass a params dict (e.g. {{'kind': 'SN', 'distance_pc': ...}}) "
                f"instead of a bare name if you want it built on demand."
            )
        raise TypeError(f"Unsupported template spec: {template!r}")

    def add_event(self, start_time_kyr, template, weight=1.0, label=None):
        """
        Add one more flux event. `template` may be an already-loaded
        FluxTemplate, the name of an already-cached one, or a params dict
        (e.g. `{"kind": "SN", "distance_pc": 200}`) -- in the last case
        the template is built via MCEq on first use and simply loaded
        from disk every time after that. Multiple events (even the same
        template at the same or nearby times, e.g. two supernovae close
        together) are superposed automatically.
        """
        resolved = self._resolve(template)
        self.events.append({
            "template": resolved, "start_time_kyr": float(start_time_kyr),
            "weight": weight, "label": label,
        })
        self._name = None  # a hand-set name can no longer describe the new content
        return resolved

    def flux(self, species, energy_gev, t_kyr):
        """
        Total flux (baseline + all active events) for one species, at
        absolute time(s) `t_kyr` (0 = present, negative = past -- see
        module docstring).
        """

        total = self.baseline.flux(species, energy_gev, t_kyr)
        total = np.broadcast_to(total, np.broadcast_shapes(np.shape(t_kyr), np.shape(energy_gev))).copy()
        for event in self.events:
            t_since = np.asarray(t_kyr, dtype=float) - event["start_time_kyr"]
            t_since *= 1.e3
            total = total + event["weight"] * event["template"].flux(species, energy_gev, t_since_kyr=t_since)
        return total

    def get_map(self, species, t_grid_kyr, e_grid_gev=None):
        """
        Returns `(t_grid_kyr, e_grid_gev, flux_map)`, `flux_map.shape ==
        (len(t_grid_kyr), len(e_grid_gev))`, ready for
        `plt.pcolormesh(e_grid, t_grid, flux_map.T)`.
        """
        if e_grid_gev is None:
            e_grid_gev = self.baseline.energy_gev
        t_grid_kyr = np.asarray(t_grid_kyr, dtype=float)
        e_grid_gev = np.asarray(e_grid_gev, dtype=float)
        tt, ee = np.meshgrid(t_grid_kyr, e_grid_gev, indexing='ij')
        return t_grid_kyr, e_grid_gev, self.flux(species, ee, tt)

try:
    from crflux.models import PrimaryFlux
except ImportError:
    PrimaryFlux = object


class SNHG12(PrimaryFlux):

    def __init__(self, model):
        self.name = 'Mod'
        self.sname = model
        self.t = model[0] * 3.154e7
        self.d = model[1] * 3.086e18
        self.params = {}
        self.SNfrac = {}
        self.rid_cutoff = {}

        mass_comp = [14, 402, 1206, 2814, 5426]
        for mcomp in mass_comp:
            self.params[mcomp] = {}
            self.SNfrac[mcomp] = 0.0

        self.rid_cutoff[1] = 4e6
        self.rid_cutoff[2] = 30e6
        self.rid_cutoff[3] = 2e9

        self.SNfrac[14] = 0.90    # H
        self.SNfrac[402] = 0.08   # He
        self.SNfrac[1206] = 0.01  # CNO
        self.SNfrac[2814] = 0.005  # MgAlSi
        self.SNfrac[5426] = 0.005  # Fe

        self.params[14][1] = (1. * 7860, 1.66, 1)    # H
        self.params[402][1] = (1. * 3550, 1.58, 2)   # He
        self.params[1206][1] = (1. * 2200, 1.63, 6)  # CNO
        self.params[2814][1] = (1. * 1430, 1.67, 14)  # MgAlSi
        self.params[5426][1] = (1. * 2120, 1.63, 26)  # Fe

        self.params[14][2] = (1 * 20, 1.4, 1)      # H
        self.params[402][2] = (1 * 20, 1.4, 2)     # He
        self.params[1206][2] = (1 * 13.4, 1.4, 6)  # CNO
        self.params[2814][2] = (1 * 13.4, 1.4, 14)  # MgAlSi
        self.params[5426][2] = (1 * 13.4, 1.4, 26)  # Fe

        self.params[14][3] = (1 * 1.7, 1.4, 1)      # H
        self.params[402][3] = (1 * 1.7, 1.4, 2)     # He
        self.params[1206][3] = (1 * 1.14, 1.4, 6)   # CNO
        self.params[2814][3] = (1 * 1.14, 1.4, 14)  # MgAlSi
        self.params[5426][3] = (1 * 1.14, 1.4, 26)  # Fe

        self.nucleus_ids = list(self.params.keys())

    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)

        flux = 0.0

        D0 = 2e28
        Q0 = 4.e52 * self.SNfrac[corsika_id]
        Ecut = 1e6
        D = D0 * (E**(1 / 3))
        Q = Q0 * (E**-2.2) * np.exp(-E / self.params[corsika_id][3][2] / Ecut)
        flux += ((Q) / (np.pi**1.5 * np.sqrt(4 * D * self.t)**3) *
                 np.exp(-self.d**2 / (4 * D * self.t))) * 1e4 * 3e10 / 4 * np.pi

        return flux

class EgalHG12(PrimaryFlux):

    def __init__(self, model):
        self.name = 'Mod'
        self.sname = model
        self.t = model[0] * 3.154e7
        self.egfactor = model[1]
        self.egcut = model[2]
        self.tcut = model[3] * 3.154e7
        self.params = {}
        self.SNfrac = {}
        self.rid_cutoff = {}

        mass_comp = [14, 402, 1206, 2814, 5426]
        for mcomp in mass_comp:
            self.params[mcomp] = {}

        self.rid_cutoff[1] = 4e6
        self.rid_cutoff[2] = 30e6
        self.rid_cutoff[3] = 2e9

        self.params[14][1] = (1. * 7860, 1.66, 1)    # H
        self.params[402][1] = (1. * 3550, 1.58, 2)   # He
        self.params[1206][1] = (1. * 2200, 1.63, 6)  # CNO
        self.params[2814][1] = (1. * 1430, 1.67, 14)  # MgAlSi
        self.params[5426][1] = (1. * 2120, 1.63, 26)  # Fe

        self.params[14][2] = (1 * 20, 1.4, 1)      # H
        self.params[402][2] = (1 * 20, 1.4, 2)     # He
        self.params[1206][2] = (1 * 13.4, 1.4, 6)  # CNO
        self.params[2814][2] = (1 * 13.4, 1.4, 14)  # MgAlSi
        self.params[5426][2] = (1 * 13.4, 1.4, 26)  # Fe

        self.params[14][3] = (1 * 1.7, 1.4, 1)      # H
        self.params[402][3] = (1 * 1.7, 1.4, 2)     # He
        self.params[1206][3] = (1 * 1.14, 1.4, 6)   # CNO
        self.params[2814][3] = (1 * 1.14, 1.4, 14)  # MgAlSi
        self.params[5426][3] = (1 * 1.14, 1.4, 26)  # Fe

        self.nucleus_ids = list(self.params.keys())

    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)

        flux = 0.0

        tdecay = np.exp(-self.t / self.tcut) if self.t > self.tcut else 1.0

        flux += self.egfactor * self.params[corsika_id][3][0] * E ** (-self.params[corsika_id][3][1] - 1.0) * \
            np.exp(-E / self.params[corsika_id][3][2] / self.egcut) * tdecay

        return flux

import numpy as np

class EnhancedHG12(PrimaryFlux):

    MODERN_MAX_CUTOFF = 15.0        
    MAX_UNSHIELDED_MULTIPLIER = 2.8

    def __init__(self, model):
        super().__init__()
        
        self.name = 'Mod'
        self.sname = model
        self.t = model[0] * 3.154e7
        self.tcut = model[2] * 3.154e7
        
        self.field_reduction_factor = float(model[1]) 

        self.params = {}
        self.rid_cutoff = {1: 4e6, 2: 30e6, 3: 2e9}
        self.z_map = {14: 1, 402: 2, 1206: 6, 2814: 14, 5426: 26}

        mass_comp = [14, 402, 1206, 2814, 5426]
        for mcomp in mass_comp:
            self.params[mcomp] = {}

        # Population 1
        self.params[14][1] = (1. * 7860, 1.66, 1)    # H
        self.params[402][1] = (1. * 3550, 1.58, 2)   # He
        self.params[1206][1] = (1. * 2200, 1.63, 6)  # CNO
        self.params[2814][1] = (1. * 1430, 1.67, 14) # MgAlSi
        self.params[5426][1] = (1. * 2120, 1.63, 26) # Fe

        # Population 2
        self.params[14][2] = (1 * 20, 1.4, 1)      
        self.params[402][2] = (1 * 20, 1.4, 2)     
        self.params[1206][2] = (1 * 13.4, 1.4, 6)  
        self.params[2814][2] = (1 * 13.4, 1.4, 14) 
        self.params[5426][2] = (1 * 13.4, 1.4, 26) 

        # Population 3
        self.params[14][3] = (1 * 1.7, 1.4, 1)      
        self.params[402][3] = (1 * 1.7, 1.4, 2)     
        self.params[1206][3] = (1 * 1.14, 1.4, 6)   
        self.params[2814][3] = (1 * 1.14, 1.4, 14)  
        self.params[5426][3] = (1 * 1.14, 1.4, 26)  

        self.nucleus_ids = list(self.params.keys())

    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)
        Z = self.z_map.get(corsika_id, 1)

        # Base HG12 flux calculation
        flux = 0.0
        for i in range(1, 4):
            p = self.params[corsika_id][i]
            flux += p[0] * E ** (-p[1] - 1.0) * np.exp(-E / (p[2] * self.rid_cutoff[i]))

        if self.t < self.tcut and self.field_reduction_factor > 1.0:
            
            turn_off_rigidity = self.MODERN_MAX_CUTOFF
            
            max_enhancement = 1.0 + (self.MAX_UNSHIELDED_MULTIPLIER - 1.0) * (1.0 - 1.0 / self.field_reduction_factor)
            
            steepness = 3.2 / np.log10(self.field_reduction_factor)
            
            E_turn = turn_off_rigidity * Z
            
            enhancement_factor = 1.0 + (max_enhancement - 1.0) / (
                1.0 + np.exp(steepness * (np.log10(E) - np.log10(E_turn)))
            )
            flux_e = flux * enhancement_factor
        else:
            flux_e = flux

        return flux_e - flux
    
def _build_template(name=None, model=None, params=None, t_since_kyr=None, interaction_model="SIBYLL2.3c"):

    import crflux.models as pm
    from MCEq.core import MCEqRun

    if model is None:
        model = pm.HillasGaisser2012
        params = "H3a"
        name = "Baseline"
        t_since_kyr = None

    if name is None:
        name = f"{model[0].__name__}_{model[1]}"

    mceq_run = MCEqRun(
    interaction_model=interaction_model,
    primary_model=(pm.HillasGaisser2012, "H3a"),
    theta_deg=0.0
    )

    e_grid = mceq_run.e_grid
    angles = DEFAULT_ANGLES_DEG

    if params is None:
        params = ()

    if t_since_kyr is None:

        mceq_run.set_primary_model(model, params)

        muons = np.zeros_like(e_grid)
        neutrons = np.zeros_like(e_grid)

        for theta in angles:
            mceq_run.set_theta_deg(theta)
            mceq_run.solve()
            muons += (mceq_run.get_solution('total_mu+', mag=0)
                    + mceq_run.get_solution('total_mu-', mag=0)) * 1e4
            neutrons += mceq_run.get_solution('n0', mag=0) * 1e4
        neutrons /= len(angles)
        muons /= len(angles)
        primary = model(params).total_flux(e_grid) * 1e4

        E0 = e_grid[32]

        g_pre = 2.6 

        anflux = (neutrons[32] / (np.exp(-1 / E0))) * np.exp(-1. / e_grid) * (e_grid / E0) ** (-g_pre)
        neutrons = np.where(e_grid > e_grid[15], neutrons, anflux)


        np.clip(neutrons, 1e-44, None, out=neutrons)
        np.clip(muons, 1e-44, None, out=muons)

    else:
        muons, neutrons, primary = [], [], []

        for age_kyr in t_since_kyr:

            age_params = (age_kyr, ) + params

            mceq_run.set_primary_model(model, age_params)

            muons_t = np.zeros_like(e_grid)
            neutrons_t = np.zeros_like(e_grid)

            for theta in angles:
                mceq_run.set_theta_deg(theta)
                mceq_run.solve()
                muons_t += (mceq_run.get_solution('total_mu+', mag=0)
                                + mceq_run.get_solution('total_mu-', mag=0))
                neutrons_t += mceq_run.get_solution('n0', mag=0)
            neutrons_t /= len(angles)
            muons_t /= len(angles)
            primary_t = model(age_params).total_flux(e_grid)

            E0 = e_grid[32]

            g_pre = 2.6 

            anflux = (neutrons_t[32] / (np.exp(-1 / E0))) * np.exp(-1. / e_grid) * (e_grid / E0) ** (-g_pre)
            neutrons_t = np.where(e_grid > e_grid[15], neutrons_t, anflux)

            np.clip(neutrons_t, 1e-44, None, out=neutrons_t)
            np.clip(muons_t, 1e-44, None, out=muons_t)

            neutrons.append(neutrons_t * 1e4)
            muons.append(muons_t * 1e4)
            primary.append(primary_t * 1e4)

    return FluxTemplate(
        name, e_grid,
        {
            "primary": np.array(primary),
            "mu+": np.array(muons) / 2.0,
            "mu-": np.array(muons) / 2.0,
            "neutron": np.array(neutrons),
        },
        t_since_kyr=t_since_kyr,
    )