"""
1-d thermal modeling functions
"""
import math

import numba as nb
import numpy as np

from mrm.shared.constants import GM, SIGMA_SB
from mrm.modeling.heat1d import planets
import mrm.modeling.heat1d.numbafied as nbf


class Model(object):

    # Initialization
    def __init__(self, planet=planets.Moon, **kwargs):
        # Initialize
        superfluous = set(kwargs).difference(set(self.__annotations__))
        if len(superfluous) > 0:
            raise TypeError(
                f"the following kwargs are not valid attributes for this "
                f"class: {superfluous}"
            )
        self._init_from_planet(planet)
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.cpcoeff = np.array(self.cpcoeff)
        self.z, self.temp = self._spatial_grid()
        self.kc = self._assign_kc()
        self._initialize_layers()
        self.sabs = self.s0 * (1.0 - self.albedo)
        self._assign_useful_parameters()
        # Model run times
        # Equilibration time -- TODO, maybe: change to convergence check
        self.equiltime = (
            self.nyearseq * self.year
            - (self.nyearseq * self.year) % self.day
        )
        # Run time for output
        self.endtime = self.equiltime + self.ndays * planet.day
        self.t = 0.0
        self.dt = self._get_time_step()
        # Check for maximum time step
        self.dtout = self.dt
        dtmax = self.day / self.nperday
        if self.dt > dtmax:
            self.dtout = dtmax
        # Array for output temperatures and local times
        self.n_steps = np.int64((self.ndays * self.day) / self.dtout)
        self.n_z = np.size(self.z)
        self.output_temp = np.zeros([self.n_steps, self.n_z])
        self.lt = np.zeros([self.n_steps])

    def _assign_useful_parameters(self):
        # useful parameter for orbit parameters
        self.useful_parameter = self.rau * (1 - self.eccentricity ** 2)
        # another useful parameter for orbit parameters
        self.another_useful_parameter = math.sqrt(self.useful_parameter * GM)
        self.sin_obliq = math.sin(self.obliquity)

    def _init_from_planet(self, planet):
        planet_map = {
            'albedo': 'albedo',
            'cpcoeff': 'cpCoeff',
            'day': 'day',
            'eccentricity': 'eccentricity',
            'emissivity': 'emissivity',
            'h': 'H',
            'ks': 'ks',
            'kd': 'kd',
            'lp': 'Lp',
            'obliquity': 'obliquity',
            'qb': 'Qb',
            'rau': 'rAU',
            'rhos': 'rhos',
            'rhod': 'rhod',
            's0': 'S',
            'year': 'year'
        }
        for attr, ref in planet_map.items():
            setattr(self, attr, getattr(planet, ref))
        self.albedo_coef_1 = planet.albedoCoef[0]
        self.albedo_coef_2 = planet.albedoCoef[1]

    def _assign_kc(self):
        ks = self.ks
        kd = self.kd
        return kd - (kd - ks) * np.exp(-self.z / self.h)

    def _initialize_layers(self):
        """initialize derived values for layers"""
        self.r350 = r350(self.chi)
        self.dz = np.diff(self.z)
        d3z = self.dz[1:] * self.dz[0:-1] * (self.dz[1:] + self.dz[0:-1])
        self.g1 = 2 * self.dz[1:] / d3z[0:]  # A.K.A. "p" in the Appendix
        self.g2 = 2 * self.dz[0:-1] / d3z[0:]  # A.K.A. "q" in the Appendix
        # this calculation differs in form from the C but gets similar results
        # -- however, it does not use the TI parameter (here called Gamma?) and
        # is most similar at TI == TI0 (55 default)
        self.rho = (
            self.rhod
            - (self.rhod - self.rhos)
            * np.exp(-self.z / self.h)
        )
        # initial conductivity profile
        self.k = nbf.thermcond(self.temp, self.kc, self.r350)
        # initial heat capacity profile
        self.cp = nbf.heat_capacity(self.cpcoeff, self.temp)

    def run(self):
        # Equilibrate the model
        # print("----equilibriating----")
        self._advance_to_equilibrium()
        # Run through end of model and store output
        self.dt = self.dtout
        self.t = 0.0  # reset simulation time
        # print("----concluding model----")
        for i in range(0, self.n_steps):
            self._advance()
            self.output_temp[i, :] = self.temp  # temperature [K]
            self.lt[i] = self.t / self.day * 24.0  # local time [hr]

    def _advance_args(self):
        return (
            self.albedo,
            self.albedo_coef_1,
            self.albedo_coef_2,
            self.another_useful_parameter,
            self.cp,
            self.cpcoeff,
            self.day,
            self.dt,
            self.dtsurf,
            self.dz,
            self.eccentricity,
            self.emissivity,
            self.equiltime,
            self.g1,
            self.g2,
            self.k,
            self.kc,
            self.lat,
            self.lp,
            self.nu,
            self.qb,
            self.r350,
            self.rau,
            self.rho,
            self.sabs,
            self.sigma,
            self.sin_obliq,
            self.t,
            self.temp,
            self.useful_parameter
        )

    def _advance_to_equilibrium(self):
        (
            self.temp,
            self.k,
            self.cp,
            self.nu,
            self.t
        ) = nbf.advance_to_equilibrium(*self._advance_args())

    def _advance(self):
        (
            self.temp,
            self.k,
            self.cp,
            self.nu
        ) = nbf.advance(*self._advance_args())
        self.t += self.dt

    def _get_time_step(self):
        dt_max = np.min(
            self.fourier_mesh_number
            * self.rho[:-1] * self.cp[:-1] * self.dz ** 2 / self.k[:-1]
        )
        return dt_max

    def _spatial_grid(self):
        """Calculate the spatial grid.

        The spatial grid is non-uniform, with layer thickness increasing
        downward.

        Returns
        -------
        np.array(float)
            Spatial node coordinates in meters.
        """
        # m: Number of layers in upper skin depth [default: 10, set in Config]
        # n: Layer increase with depth:
        #   dz[i] = dz[i-1]*(1+1/n) [default: 5, set in Config]
        # b: Number of skin depths to bottom layer [default: 20, set in Config]
        # ztscale: arbitrary (?) scaling parameter for temperature rough-in

        # starting temperature values
        # surface temperature
        t0s = t_radeq(self.albedo, self.s0, self.lat, self.rau)
        t0d = t0s / (2 ** 0.5)  # subsurface temperature
        # NOTE: adding temperature-dependent skin depth from C code
        # NOTE: planet.ks is similar to kc[0] in C code
        zs = skin_depth(self.cpcoeff, self.ks, self.rhos, self.day, t0s)
        dz = np.zeros(1) + zs / self.nskin  # thickness of uppermost layer
        z = np.zeros(1)
        zmax = zs * self.nskinbot  # depth of deepest layer
        # initialize geometrically-increasing layer depth / thickness
        i = 0
        while z[i] < zmax:
            i += 1
            h = dz[i - 1] * (1 + 1 / self.layer_scale_factor)
            dz = np.append(dz, h)
            # NOTE: corrected apparently-accidental dz[i] reference
            z = np.append(z, z[i - 1] + dz[i - 1])  # depth of layer i
        # TODO: in the C code, the final layer is everywhere ignored, used only
        #  for thickness (the bottom layer for which temperatures are computed
        #  is nlayers - 1) -- we may clean this up for clarity later
        temp = []
        # temperature curve rough-in, equvalent to middle section of tiProfile
        # in C version -- necessary to make model stable at certain settings
        # (and also speeds it up)
        for depth in z:
            temp.append(t0d - (t0d - t0s) * np.exp(-depth / self.ztscale))
        return z, np.array(temp)

    albedo: float
    albedo_coef_1: float
    albedo_coef_2: float
    another_useful_parameter: float
    # Radiative conductivity parameter [Mitchell and de Pater, 1994]
    chi: float = 2.7
    cp: nb.float64[:]
    cpcoeff: nb.float64[:]
    day: float
    dec = np.float64()  # solar declination [rad]
    dt: float
    # surface temperature accuracy (proportional to absolute temperature)
    dtsurf: float = 0.01
    dz: nb.float64[:]
    eccentricity: float
    emissivity: float
    equiltime: float
    fourier_mesh_number: float = 0.49  # must be <= 0.5 for stability
    g1: nb.float64[:]
    g2: nb.float64[:]
    h: np.float64()
    k: nb.float64[:]
    kc: nb.float64[:]
    ks: np.float64()  # Solid (phonon) conductivity at surface [W.m-1.K-1]
    kd: np.float64()  # Solid (phonon) conductivity at depth [W.m-1.K-1]
    lat: float = 0
    lp: float
    # Layer increase with depth: dz[i] = dz[i-1]*(1+1/layer_scale_factor)
    #  -- elsewhere NN, n
    layer_scale_factor: int = 10
    ndays: int = 1  # number of days to ouptut
    nperday: int = 24  # minimum time steps per diurnal cycle
    nskin: int = 10  # Number of layers in upper skin depth, elsewhere 'm'
    nskinbot: int = 20  # Number of skin depths to bottom layer, elsewhere 'b'
    nyearseq: int = 1  # equilibration time [orbits]
    nu: float = np.float64()  # orbital true anomaly [rad]
    nudot = np.float64()  # rate of change of true anomaly [rad/s]
    obliquity: np.float64()
    qb: float
    r350: float
    rau: float  # solar distance [AU]
    rho: nb.float64[:]
    rhod: np.float64()
    rhos: np.float64()
    s0: float  # Solar constant at 1 AU [W.m-2]
    sabs: float
    sigma: float = SIGMA_SB  # Stefan-Boltzmann Constant, 5.67051196e-8
    sin_obliq: float
    temp: nb.float64[:]
    t: float
    useful_parameter: float
    year: float
    ztscale: float = 0.1  # parameter for temperature initialization scaling


def t_radeq(albedo, s0, lat, rau):
    """Calculate radiative equilibrium temperature at local noontime."""
    # NOTE: modified to add distance term
    # NOTE: retains addition of an emissivity term absent from equivalent
    # initializing operation in C, but it makes differences in results
    # only at the thousandth-of-a-degree level for the Moon
    return (
        (1 - albedo) / (SIGMA_SB * rau ** 2) * s0 * np.cos(lat)
    ) ** 0.25


def skin_depth(cpcoeff, ks, rhos, day, t0s):
    """Calculate Thermal skin depth [m]."""
    # NOTE: adding temperature-dependent skin depth from C code
    cp0 = nbf.heat_capacity_scalar(cpcoeff, t0s)
    # thermal diffusivity = k / (rho * cp)[m2.s - 1]
    # NOTE: ks is not treated as a constant in the C code, but as dependent on
    # TI. However, at the default value of 55, it is very very close to 0.0008.
    kappa = ks / (rhos * cp0)
    return np.sqrt(kappa * day / np.pi)


def r350(chi):
    return chi / 350 ** 3
