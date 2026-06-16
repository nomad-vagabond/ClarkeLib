"""
ClarkeLib - an open library for calculation of basic space elevator
parameters and loads.
With respect and gratitude to Sir Arthur C. Clarke.

Copyright 2015 Vadym Pasko all rights reserved,
Vadym Pasko <keenon3d@gmail.com>
Permission to use, modify, and distribute this software is given under the
terms of the GPL 2 license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.

Units of all physical parameters are in SI (International System of Units),
except Specific Strength of the tether, which uses MYuri - widely spread unit
in literature on space elevator.

Numerical Python port of clarkelib.sage.
Symbolic SageMath operations replaced with scipy.integrate.quad and
scipy.optimize (brentq, fsolve).
"""

import math
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq, fsolve

# ---------------------------------------------------------------------------
# Error / warning messages (identical to clarkelib.sage)
# ---------------------------------------------------------------------------

_warn_msg = {
    1: "since you've specified hcw, khcw parameter is ignored"
}

_error_msg = {
    0: "'planet' parameter must be string.",
    1: "'anchor_safe' parameter must be 'heavy' or 'light'.",
    2: "taper_func parameter must be a list or tuple of tapering fuctions names.",
    3: "Integral has not been found. Try use another tapering function.",
    4: ("No solution was found. "
        "Probably material is too weak for the selected planet."
        " \n Try using better material."),
    5: ("'material' parameter must be a dictionary "
        "with 'rho', 'ss' and 'ksafe' keys with numeric values."),
    6: ("Initial value of the counterweight altitude must "
        "be > then altitude of the planet's synchronous orbit."),
    7: ("Final value of the counterweight altitude must "
        "be > then altitude of the planet's synchronous orbit."),
    8: "hcw must be larger then altitude of the planet's synchronous orbit",
    9: "'planet' parameter must be one of the available planet names: ",
}

# ---------------------------------------------------------------------------
# Taper-function factories
#
# Each factory returns a callable  f(h) -> dimensionless ratio A(h)/A_anchor.
# Low part  (R  <= h <= hsyn): f(R)=1,    f(hsyn)=ksyn
# High part (hsyn <= h <= hcw): f(hsyn)=ksyn, f(hcw)=kcw
# ---------------------------------------------------------------------------

def _make_exp_taper(GM, hsyn, R, rho, ksafe, ts):
    """Exponential (constant-stress) taper, normalised so that f(R) = 1."""
    C = rho * ksafe / ts
    ref = GM * C * (1.0/R + R**2 / (2.0*hsyn**3))

    return lambda h: math.exp(ref - GM * C * (1.0/h + h**2 / (2.0*hsyn**3)))


def _make_linear_low(ksyn, R, hsyn):
    return lambda h: (hsyn - h) / (hsyn - R) + ksyn * (h - R) / (hsyn - R)


def _make_sqrt_low(ksyn, R, hsyn):
    return lambda h: (ksyn - 1) * math.sqrt(h - R) / math.sqrt(hsyn - R) + 1.0


def _make_cbrt_low(ksyn, R, hsyn):
    return lambda h: (ksyn - 1) * (h - R)**(1.0/3) / (hsyn - R)**(1.0/3) + 1.0


def _make_quadratic_low(ksyn, R, hsyn):
    return lambda h: -(ksyn - 1) * (h - hsyn)**2 / (R - hsyn)**2 + ksyn


def _make_cubic_low(ksyn, R, hsyn):
    return lambda h: -(ksyn - 1) * (h - hsyn)**3 / (R - hsyn)**3 + ksyn


def _make_quartic_low(ksyn, R, hsyn):
    return lambda h: -(ksyn - 1) * (h - hsyn)**4 / (R - hsyn)**4 + ksyn


def _make_linear_high(ksyn, kcw, hsyn, hcw):
    return lambda h: kcw * (h - hsyn) / (hcw - hsyn) + ksyn * (hcw - h) / (hcw - hsyn)


def _make_quadratic_high(ksyn, kcw, hsyn, hcw):
    return lambda h: -(ksyn - kcw) * (h - hsyn)**2 / (hcw - hsyn)**2 + ksyn


def _make_cubic_high(ksyn, kcw, hsyn, hcw):
    return lambda h: -(ksyn - kcw) * (h - hsyn)**3 / (hcw - hsyn)**3 + ksyn


def _make_quartic_high(ksyn, kcw, hsyn, hcw):
    return lambda h: -(ksyn - kcw) * (h - hsyn)**4 / (hcw - hsyn)**4 + ksyn


_taperfunc_low_makers = {
    'linear':    _make_linear_low,
    'sqrt':      _make_sqrt_low,
    'cbrt':      _make_cbrt_low,
    'quadratic': _make_quadratic_low,
    'cubic':     _make_cubic_low,
    'quartic':   _make_quartic_low,
}

_taperfunc_high_makers = {
    'linear':    _make_linear_high,
    'quadratic': _make_quadratic_high,
    'cubic':     _make_cubic_high,
    'quartic':   _make_quartic_high,
}

_taperfunc_low_names  = list(_taperfunc_low_makers.keys())  + ['exp']
_taperfunc_high_names = list(_taperfunc_high_makers.keys()) + ['exp']


# ---------------------------------------------------------------------------
# Planet
# ---------------------------------------------------------------------------

class Planet(object):
    """
    Planet of the solar system. Currently implemented Mars and Earth.

    Parameters
    ----------
    name : str
           Name pointing to one of available planets.
    """
    G     = float(6.67384e-11)
    names = ("Earth", "Mars")

    def __init__(self, name):
        if type(name) is not str:
            raise ValueError(_error_msg[0])
        if name not in self.names:
            planetlist = "'%s', " * (len(self.names) - 1) + "'%s'."
            raise ValueError(_error_msg[9] + planetlist % self.names)
        self.name = name
        if name == "Earth":
            self._setEarth()
        elif name == "Mars":
            self._setMars()
        self._setSynOrb()

    def _setEarth(self):
        self.M     = float(5.97219e24)
        self.R     = float(6371e3)
        self.OMEGA = float(7.292e-5)

    def _setMars(self):
        self.M     = float(6.4185e23)
        self.R     = float(3396e3)
        self.OMEGA = float(7.077651712e-5)

    def _setSynOrb(self):
        self.GM   = self.G * self.M
        self.HSYN = (self.GM / self.OMEGA**2) ** (1.0 / 3)


# ---------------------------------------------------------------------------
# Counterweight
# ---------------------------------------------------------------------------

class Counterweight(object):
    """Space Elevator counterweight."""

    def __init__(self, hcw, planet):
        self.planet   = planet
        self.hcw      = hcw
        self.mcw      = None
        self.force_cw = None

    def getCounterweightForce(self, symbolic=False):
        """Compute net upward centrifugal force minus gravity for the counterweight."""
        pl = self.planet
        self.force_cw = (self.mcw * pl.OMEGA**2 * self.hcw
                         - self.mcw * pl.GM / self.hcw**2)


# ---------------------------------------------------------------------------
# Tether
# ---------------------------------------------------------------------------

class Tether(object):
    """
    Space Elevator tether.
    Attributes with _low suffix: lower part (R <= h <= HSYN).
    Attributes with _high suffix: upper part (HSYN <= h <= hcw).

    Parameters
    ----------
    rho          : float  – tether material density, kg/m³
    ss           : float  – specific strength, MYuri
    ksafe        : float  – safety factor
    planet       : Planet
    cw           : Counterweight
    pull_force   : float  – anchor_force + total climber weight, N
    anchor_force : float  – force at anchor point, N
    taper_func   : tuple  – (low_name, high_name)
    climbers     : list   – [(force, h_start, h_end), ...]
    anchor_safe  : str    – 'heavy' or 'light'
    """

    def __init__(self, rho, ss, ksafe, planet, cw,
                 pull_force, anchor_force, taper_func, climbers, anchor_safe):

        self.rho        = float(rho)
        self.ss         = float(ss)
        self.ts         = self.ss * rho * 1e6       # tensile strength, Pa
        self.ksafe      = float(ksafe)
        self.sigma_work = self.ts / self.ksafe       # working stress, Pa

        self.planet      = planet
        self.cw          = cw
        self.pull_force  = pull_force
        self.anchor_force = anchor_force
        self.climbers    = climbers
        self.anchor_safe = anchor_safe

        if anchor_safe == 'heavy':
            anchor_pullforce = pull_force
        elif anchor_safe == 'light':
            anchor_pullforce = anchor_force + climbers[0][0]
        else:
            raise ValueError(_error_msg[1])

        # f0 [kg/m]: linear mass density at anchor point
        # A_anchor = f0 / rho [m²]: cross-section area at anchor
        self.f0       = self.rho * self.ksafe * anchor_pullforce / self.ts
        self.A_anchor = self.f0 / self.rho

        if type(taper_func) in (list, tuple):
            self.tf_low  = taper_func[0]
            self.tf_high = taper_func[1]
        else:
            raise ValueError(_error_msg[2])

        self.solveTether()

    # ------------------------------------------------------------------
    # Taper-function builders
    # ------------------------------------------------------------------

    def _setTaperFunctions(self):
        """Validate taper function names."""
        if self.tf_low not in _taperfunc_low_names:
            raise ValueError("tapering function %s has not been found. "
                             "Try one of: %s" % (self.tf_low,
                             ', '.join(_taperfunc_low_names)))
        if self.tf_high not in _taperfunc_high_names:
            raise ValueError("Tapering function %s has not been found. "
                             "Try one of: %s" % (self.tf_high,
                             ', '.join(_taperfunc_high_names)))

    def _makeExpTaper(self):
        pl = self.planet
        return _make_exp_taper(pl.GM, pl.HSYN, pl.R,
                               self.rho, self.ksafe, self.ts)

    def _buildTaperLow(self, ksyn):
        """Return dimensionless taper callable for the lower part."""
        pl = self.planet
        if self.tf_low == 'exp':
            return self._makeExpTaper()
        return _taperfunc_low_makers[self.tf_low](ksyn, pl.R, pl.HSYN)

    def _buildTaperHighExp(self, ksyn_fix):
        """Return exp high taper rescaled by ksyn_fix for continuity at HSYN."""
        exp_fn = self._makeExpTaper()
        def taper_high(h, _e=exp_fn, _f=ksyn_fix):
            return _e(h) * _f
        return taper_high

    def _buildTaperHighPoly(self, ksyn, kcw):
        """Return polynomial high taper callable."""
        pl  = self.planet
        hcw = self.cw.hcw
        return _taperfunc_high_makers[self.tf_high](ksyn, kcw, pl.HSYN, hcw)

    # ------------------------------------------------------------------
    # Force integration (all force_* attributes in m²/s² — dimensionless
    # taper × net gravity integrated over altitude.
    # Multiply by A_anchor * rho to get Newtons.)
    # ------------------------------------------------------------------

    def _findForcesLow(self, taper_fn):
        pl = self.planet
        def grav_int(h): return taper_fn(h) * pl.GM / h**2
        def cent_int(h): return taper_fn(h) * pl.OMEGA**2 * h
        self.gravforce_low, _ = quad(grav_int, pl.R, pl.HSYN, limit=100)
        self.centforce_low, _ = quad(cent_int, pl.R, pl.HSYN, limit=100)
        self.force_low = self.gravforce_low - self.centforce_low

    def _findForcesHigh(self, taper_fn):
        pl  = self.planet
        hcw = self.cw.hcw
        def grav_int(h): return taper_fn(h) * pl.GM / h**2
        def cent_int(h): return taper_fn(h) * pl.OMEGA**2 * h
        self.gravforce_high, _ = quad(grav_int, pl.HSYN, hcw, limit=100)
        self.centforce_high, _ = quad(cent_int, pl.HSYN, hcw, limit=100)
        self.force_high  = self.gravforce_high - self.centforce_high
        self.gravforce   = self.gravforce_low + self.gravforce_high
        self.centforce   = self.centforce_low + self.centforce_high

    # ------------------------------------------------------------------
    # Solver methods
    # ------------------------------------------------------------------

    def _getTetherAreas(self, ksyn, kcw, part='both'):
        A0 = self.A_anchor
        if part in ('both', 'low'):
            self.s1 = A0               # anchor cross-section, m²
            self.s2 = A0 * ksyn        # GEO cross-section, m²
        if part in ('both', 'high'):
            self.s3 = A0 * kcw         # CW cross-section, m²

    def _findKsyn(self, exp_fn):
        """
        Compute self.ksyn (and self.ksyn_fix for non-exp/exp case).

        For exp lower taper: ksyn = exp_fn(HSYN), ksyn_fix = 1.
        For polynomial lower taper: solve stress-balance at GEO with brentq.
        """
        pl = self.planet
        if self.tf_low == 'exp':
            self.ksyn     = exp_fn(pl.HSYN)
            self.ksyn_fix = 1.0
            return

        A0 = self.A_anchor

        def residual(ksyn_v):
            tl = self._buildTaperLow(ksyn_v)
            self._findForcesLow(tl)
            F_geo = self.pull_force + self.force_low * A0 * self.rho
            return F_geo / (A0 * ksyn_v) - self.sigma_work

        self.ksyn = brentq(residual, 0.01, 500.0, xtol=1e-8, maxiter=200)

        if self.tf_high == 'exp':
            self.ksyn_fix = self.ksyn / exp_fn(pl.HSYN)
        else:
            self.ksyn_fix = 1.0

    def _findKcwMcw(self, taper_low, taper_high_builder):
        """
        Solve the 2-equation system for kcw and mcw.

        taper_high_builder(kcw_v) returns a taper callable for the high part.
        This allows the integrand to depend on kcw for polynomial high tapers.

        Equation 1 (equilibrium):
            force_low_N + force_high_N + pull_force = m_cw * cw_spec
        Equation 2 (stress at CW altitude):
            m_cw * cw_spec / (A_anchor * kcw) = sigma_work
        """
        pl  = self.planet
        A0  = self.A_anchor

        # Low-part forces are fixed; compute once before the solve loop.
        self._findForcesLow(taper_low)
        force_low_N = self.force_low * A0 * self.rho

        cw_spec = pl.OMEGA**2 * self.cw.hcw - pl.GM / self.cw.hcw**2

        def system(x):
            kcw_v, mcw_v = float(x[0]), float(x[1])
            th = taper_high_builder(kcw_v)
            self._findForcesHigh(th)
            force_high_N = self.force_high * A0 * self.rho
            F_cw = mcw_v * cw_spec
            s3   = A0 * kcw_v
            eq1  = force_low_N + force_high_N + self.pull_force - F_cw
            eq2  = F_cw / s3 - self.sigma_work
            return [eq1, eq2]

        sol, _, ier, msg = fsolve(system, [1.0, 1e6], full_output=True)
        if ier != 1:
            raise ValueError(_error_msg[4])

        self.kcw      = float(sol[0])
        self.cw.mcw   = float(sol[1])

        # Ensure final force attributes match the converged solution.
        self._findForcesHigh(taper_high_builder(self.kcw))

    def solveTether(self):
        """Find tether parameters for the current counterweight altitude."""

        self._setTaperFunctions()
        exp_fn = self._makeExpTaper()

        if self.tf_low != 'exp' and self.tf_high == 'exp':
            # ---------------------------------------------------------------
            # Polynomial lower / exponential upper
            # ---------------------------------------------------------------
            # 1. Find ksyn from stress balance at GEO (uses only lower part).
            self._findKsyn(exp_fn)           # sets self.ksyn, self.ksyn_fix

            # 2. Build final lower taper with converged ksyn.
            taper_low = self._buildTaperLow(self.ksyn)
            self._findForcesLow(taper_low)
            self._getTetherAreas(self.ksyn, 1.0, part='low')

            # 3. Upper exp taper rescaled for cross-section continuity at HSYN.
            _fix = self.ksyn_fix
            def taper_high_builder(kcw_v, _f=_fix):
                return self._buildTaperHighExp(_f)

        else:
            # ---------------------------------------------------------------
            # exp/exp  |  polynomial/polynomial  |  exp/polynomial
            # ---------------------------------------------------------------
            self._findKsyn(exp_fn)           # sets self.ksyn, self.ksyn_fix
            taper_low = self._buildTaperLow(self.ksyn)
            self._findForcesLow(taper_low)

            if self.tf_high == 'exp':
                _fix = self.ksyn_fix
                def taper_high_builder(kcw_v, _f=_fix):
                    return self._buildTaperHighExp(_f)
            else:
                _ksyn = self.ksyn
                def taper_high_builder(kcw_v, _k=_ksyn):
                    return self._buildTaperHighPoly(_k, kcw_v)

        # Solve for kcw and mcw.
        self._findKcwMcw(taper_low, taper_high_builder)

        # Validate
        if self.ksyn < 0 or self.kcw < 0 or self.f0 < 0 or self.cw.mcw < 0:
            raise ValueError(_error_msg[4])

        # Store final taper callables.
        self._taper_low_fn  = taper_low
        self._taper_high_fn = taper_high_builder(self.kcw)

        self._getTetherAreas(self.ksyn, self.kcw)
        self.cw.getCounterweightForce(symbolic=False)

        # Recompute forces with final tapers.
        self._findForcesLow(self._taper_low_fn)
        self._findForcesHigh(self._taper_high_fn)
        self.tetherforce_low  = self.force_low  * self.A_anchor * self.rho
        self.tetherforce_high = self.force_high * self.A_anchor * self.rho

    # ------------------------------------------------------------------
    # Post-solve helpers
    # ------------------------------------------------------------------

    def getTetherMass(self):
        """Compute tether mass by numerical integration."""
        pl  = self.planet
        A0  = self.A_anchor
        rho = self.rho
        mass_low,  _ = quad(lambda h: self._taper_low_fn(h)  * A0 * rho,
                            pl.R, pl.HSYN, limit=100)
        mass_high, _ = quad(lambda h: self._taper_high_fn(h) * A0 * rho,
                            pl.HSYN, self.cw.hcw, limit=100)
        self.m_tether = mass_low + mass_high

    def getTaperSectionLaws(self):
        """Store section-area functions [m²] for profile plotting."""
        A0 = self.A_anchor
        _tl = self._taper_low_fn
        _th = self._taper_high_fn
        self.sections_low  = lambda h: _tl(h) * A0
        self.sections_high = lambda h: _th(h) * A0

    def getFoceDistribution(self, pull_force):
        """
        Build callables fd_low(h) and fd_high(h) for inner-force distribution.
        pull_force: tension at the anchor (N).
        """
        pl  = self.planet
        A0  = self.A_anchor
        rho = self.rho
        _tl = self._taper_low_fn
        _th = self._taper_high_fn
        _tfl = self.tetherforce_low

        def _fd_low(h):
            val, _ = quad(
                lambda hp: _tl(hp) * A0 * rho * (pl.GM/hp**2 - pl.OMEGA**2*hp),
                pl.R, h, limit=100)
            return pull_force + val

        def _sd_low(h):
            return _fd_low(h) / (_tl(h) * A0)

        def _fd_high(h):
            val, _ = quad(
                lambda hp: _th(hp) * A0 * rho * (pl.GM/hp**2 - pl.OMEGA**2*hp),
                pl.HSYN, h, limit=100)
            return pull_force + _tfl + val

        def _sd_high(h):
            return _fd_high(h) / (_th(h) * A0)

        self.fd_low  = _fd_low
        self.sd_low  = _sd_low
        self.fd_high = _fd_high
        self.sd_high = _sd_high

    def getForceDistWithClimbers(self, climbers, pull_force):
        """
        Build per-climber-segment force/stress distribution functions
        for the lower tether part.
        """
        pl  = self.planet
        A0  = self.A_anchor
        rho = self.rho
        _tl = self._taper_low_fn

        self.climb_forcedist_low  = []
        self.climb_stressdist_low = []
        self.climber_sections     = []

        cumulative_pull = self.anchor_force
        for climber in climbers:
            climber_force, h_start, h_end = climber
            cumulative_pull += climber_force
            _pull = cumulative_pull  # capture by value

            def _fd(h, _p=_pull):
                val, _ = quad(
                    lambda hp: _tl(hp) * A0 * rho * (pl.GM/hp**2 - pl.OMEGA**2*hp),
                    pl.R, h, limit=100)
                return _p + val

            def _sd(h, _p=_pull):
                return _fd(h, _p) / (_tl(h) * A0)

            self.climb_forcedist_low.append(_fd)
            self.climb_stressdist_low.append(_sd)
            self.climber_sections.append((h_start, h_end))

    def getEquilibriumAccuracy(self):
        """Return equilibrium residual (should be ~0), N."""
        return (self.tetherforce_low + self.tetherforce_high
                + self.pull_force - self.cw.force_cw)


# ---------------------------------------------------------------------------
# SpaceElevator
# ---------------------------------------------------------------------------

class SpaceElevator(object):
    """
    Space Elevator main class.

    Evaluates all main Space Elevator design parameters for a given
    counterweight altitude hcw.

    Parameters
    ----------
    n_climbers   : int, optional
    m_climber    : float, optional  – single climber mass, kg
    anchor_force : float, optional  – anchor point force, N
    planet       : str, optional
    material     : dict or 'default'
    khcw         : float or None  – ratio of CW altitude (above surface)
                                    to HSYN (above surface)
    hcw          : float or None  – CW altitude above surface, m
    taper_func   : tuple  – ('low_name', 'high_name')
    anchor_safe  : str    – 'heavy' or 'light'
    """

    def __init__(self, n_climbers=1, m_climber=0, anchor_force=100,
                 planet='Earth', material='default', khcw=None, hcw=None,
                 taper_func=('exp', 'exp'), anchor_safe='heavy'):

        n_climbers   = self._checkNumeric('n_climbers', n_climbers, int)
        self.n_climbers = max(1, n_climbers)
        m_climber    = self._checkNumeric('m_climber', m_climber, float)
        self.m_climber = m_climber
        anchor_force = self._checkNumeric('anchor_force', anchor_force, float)

        if khcw is not None:
            khcw = self._checkNumeric('khcw', khcw, float)
        if hcw is not None:
            hcw = self._checkNumeric('hcw', hcw, float)

        if material == 'default':
            rho, ss, ksafe = 1300.0, 45.0, 1.3
        else:
            try:
                rho   = float(material['rho'])
                ss    = float(material['ss'])
                ksafe = float(material['ksafe'])
            except Exception:
                raise ValueError(_error_msg[5])

        self.planet = Planet(planet)
        self.getClimbersForce()
        pull_force = anchor_force + self.climbers_force_total

        if hcw is None:
            if khcw is not None:
                hcw_abs = (self.planet.HSYN - self.planet.R) * khcw + self.planet.R
            else:
                hcw_abs = (self.planet.HSYN - self.planet.R) * 2 + self.planet.R
        else:
            hcw_abs = hcw + self.planet.R
            if khcw is not None:
                warnings.warn(_warn_msg[1])
            if hcw_abs < self.planet.HSYN:
                raise ValueError(_error_msg[8])

        self.cw = Counterweight(hcw_abs, self.planet)
        self.tether = Tether(rho, ss, ksafe, self.planet, self.cw,
                             pull_force, anchor_force, taper_func,
                             self.climbers, anchor_safe)

        self.m0, self.y0 = self._getMassCenter(self.tether, self.cw)
        self.plots = ResplotSE(self)

    def _checkNumeric(self, name, value, value_type):
        erroer_msg = 'type of %s parameter must be numeric' % name
        try:
            return value_type(value)
        except Exception:
            raise ValueError(erroer_msg)

    def _getSingleClimberWeight(self, hclimber):
        return (self.m_climber * self.planet.GM / hclimber**2
                - self.m_climber * self.planet.OMEGA**2 * hclimber)

    def getClimbersForce(self):
        pl = self.planet
        hspace = np.linspace(pl.R, pl.HSYN, self.n_climbers + 1)
        self.climbers_force_total = 0.0
        self.climbers = []
        for i, hclimber in enumerate(hspace[:-1]):
            climber_force = self._getSingleClimberWeight(hclimber)
            topbound = hspace[i + 1]
            self.climbers.append((climber_force, hclimber, topbound))
            self.climbers_force_total += climber_force

    def _getMassCenter(self, th, cw):
        th.getTetherMass()
        pl = self.planet
        m0 = th.m_tether + cw.mcw

        cw_moment           = cw.mcw * (cw.hcw - pl.HSYN)
        tether_high_moment  = (cw.hcw - pl.HSYN)**2 * (th.s2 + 2*th.s3) * th.rho / 6
        tether_low_moment   = (pl.HSYN - pl.R)**2   * (th.s2 + 2*th.s1) * th.rho / 6
        se_moment = cw_moment + tether_high_moment - tether_low_moment
        y0        = se_moment / m0

        return m0, y0

    def setHcw(self, hcw, check=True):
        hcw_abs = hcw + self.planet.R
        if check:
            hcw_abs = float(hcw_abs)
        if hcw_abs < self.planet.HSYN:
            raise ValueError(_error_msg[8])
        self.cw          = Counterweight(hcw_abs, self.planet)
        self.tether.cw   = self.cw
        self.tether.solveTether()
        self.m0, self.y0 = self._getMassCenter(self.tether, self.cw)
        self.plots.update(self)

    def varyCWAltitude(self, hcw_start, hcw_end, hcw_num):
        """Find SE target parameters in a range of counterweight altitudes."""
        hcw_start = self._checkNumeric('hcw_start', hcw_start, float)
        hcw_end   = self._checkNumeric('hcw_end',   hcw_end,   float)
        hcw_num   = self._checkNumeric('hcw_num',   hcw_num,   int)

        if hcw_start < self.planet.HSYN:
            raise ValueError(_error_msg[6])
        if hcw_end < self.planet.HSYN:
            raise ValueError(_error_msg[7])

        self.plots.cleanPointLists()
        hspace = np.linspace(hcw_start, hcw_end, hcw_num + 1)[:-1]

        m0_min    = None
        hcw_m0min = None
        for hcw in hspace:
            self.setHcw(hcw, check=False)
            if m0_min is None or self.m0 < m0_min:
                m0_min    = self.m0
                hcw_m0min = self.cw.hcw
            self.plots.hcw_points.append(hcw)
            self.plots.m0_points.append(self.m0)
            self.plots.mcw_points.append(self.cw.mcw)
            self.plots.mtether_points.append(self.tether.m_tether)

        self.plots.update(self)
        self.plots.hcw_max = hcw

        return hcw_m0min - self.planet.R


# ---------------------------------------------------------------------------
# ResplotSE
# ---------------------------------------------------------------------------

class ResplotSE(object):
    """
    Collection of plots for Space Elevator parameter analyses.

    Notes
    -----
    Altitude on all plots is measured from the planet's surface.
    """

    def __init__(self, parent):
        self.cleanPointLists()
        self.update(parent)

    def update(self, parent):
        self.se     = parent
        self.planet = parent.planet
        self.tether = parent.tether
        self.cw     = parent.cw

    def cleanPointLists(self):
        self.hcw_points    = []
        self.m0_points     = []
        self.mcw_points    = []
        self.mtether_points = []

    # ------------------------------------------------------------------
    # Point-collection helpers
    # ------------------------------------------------------------------

    def _getDistPoints(self, forcedist_fn, stressdist_fn, hstart, hend, pnum=50):
        """Evaluate force [kN] and stress [GPa] callables over altitude range."""
        ap, fp, sp = [], [], []
        for h in np.linspace(hstart, hend, pnum):
            ap.append((h - self.planet.R) * 1e-6)
            fp.append(forcedist_fn(h)  * 1e-3)
            sp.append(stressdist_fn(h) * 1e-9)
        return ap, fp, sp

    def _getDistPointsExp(self, taper_fn, sections_fn, constforce,
                          hstart, hend, hzero, pnum=50):
        """
        Evaluate force/stress for an exp-type taper by integrating from hzero.
        constforce : tension at hzero, N.
        taper_fn   : dimensionless taper callable.
        sections_fn: area callable, m².
        """
        pl  = self.planet
        th  = self.tether
        A0  = th.A_anchor
        rho = th.rho

        ap, fp, sp = [], [], []
        for h in np.linspace(hstart, hend, pnum):
            grav, _ = quad(lambda hp: taper_fn(hp) * A0 * rho * pl.GM / hp**2,
                           hzero, h, limit=100)
            cent, _ = quad(lambda hp: taper_fn(hp) * A0 * rho * pl.OMEGA**2 * hp,
                           hzero, h, limit=100)
            force   = constforce + grav - cent
            stress  = force / sections_fn(h)
            ap.append((h - pl.R) * 1e-6)
            fp.append(force  * 1e-3)
            sp.append(stress * 1e-9)
        return ap, fp, sp

    def _getProfilePoints(self, section_fn, hstart, hend, pnum=50):
        ap, pp = [], []
        for h in np.linspace(hstart, hend, pnum):
            ap.append((h - self.planet.R) * 1e-6)
            pp.append(section_fn(h) * 1e6)   # m² → mm²
        return ap, pp

    # ------------------------------------------------------------------
    # Arrangement helpers
    # ------------------------------------------------------------------

    def _arrangeFSPointsLow(self, climbers):
        th = self.tether
        pl = self.planet
        se = self.se

        if th.tf_low == 'exp':
            constforce = th.anchor_force
            if climbers:
                ap_low, fp_low, sp_low = [], [], []
                for climbforce, hstart, hend in se.climbers:
                    constforce += climbforce
                    ap, fp, sp = self._getDistPointsExp(
                        th._taper_low_fn, th.sections_low,
                        constforce, hstart, hend, pl.R)
                    ap_low += ap; fp_low += fp; sp_low += sp
            else:
                ap_low, fp_low, sp_low = self._getDistPointsExp(
                    th._taper_low_fn, th.sections_low,
                    self.pullforce, pl.R, pl.HSYN, pl.R)
        else:
            if climbers:
                ap_low, fp_low, sp_low = [], [], []
                for fd, sd, (h_start, h_end) in zip(
                        th.climb_forcedist_low,
                        th.climb_stressdist_low,
                        th.climber_sections):
                    ap, fp, sp = self._getDistPoints(fd, sd, h_start, h_end)
                    ap_low += ap; fp_low += fp; sp_low += sp
            else:
                ap_low, fp_low, sp_low = self._getDistPoints(
                    th.fd_low, th.sd_low, pl.R, pl.HSYN)

        return ap_low, fp_low, sp_low

    def _arrangeFSPointsHigh(self):
        th  = self.tether
        pl  = self.planet
        cw  = self.cw

        if th.tf_high == 'exp':
            constforce = self.pullforce + th.tetherforce_low
            ap_high, fp_high, sp_high = self._getDistPointsExp(
                th._taper_high_fn, th.sections_high,
                constforce, pl.HSYN, cw.hcw, pl.HSYN)
        else:
            ap_high, fp_high, sp_high = self._getDistPoints(
                th.fd_high, th.sd_high, pl.HSYN, cw.hcw)

        return ap_high, fp_high, sp_high

    def _arrangeFSPoints(self, climbers):
        th = self.tether
        se = self.se

        th.getTaperSectionLaws()
        th.getFoceDistribution(self.pullforce)
        if th.tf_low != 'exp':
            th.getForceDistWithClimbers(se.climbers, self.pullforce)

        ap_low,  fp_low,  sp_low  = self._arrangeFSPointsLow(climbers)
        ap_high, fp_high, sp_high = self._arrangeFSPointsHigh()

        altitude_points = [ap_low[0]]  + ap_low  + ap_high + [ap_high[-1]]
        force_points    = [0.0]        + fp_low  + fp_high + [0.0]
        stress_points   = [0.0]        + sp_low  + sp_high + [0.0]

        return altitude_points, force_points, stress_points

    def _arrangeTetherProfilePoints(self):
        th = self.tether
        pl = self.planet
        cw = self.cw

        th.getTaperSectionLaws()

        ap_low,  pp_low  = self._getProfilePoints(th.sections_low,  pl.R,    pl.HSYN)
        ap_high, pp_high = self._getProfilePoints(th.sections_high, pl.HSYN, cw.hcw)

        ap = [ap_low[0]]  + ap_low  + ap_high + [ap_high[-1]]
        pp = [0.0]        + pp_low  + pp_high + [0.0]
        return ap, pp

    # ------------------------------------------------------------------
    # Public plot methods
    # ------------------------------------------------------------------

    def mplplotMvsHcw(self, figsize=(8, 10)):
        """Plot SE component masses vs counterweight altitude."""
        pl = self.planet
        plt.figure(figsize=figsize)
        x = [(p - pl.R) * 1e-6 for p in self.hcw_points]

        plt.subplot(3, 1, 1)
        plt.plot(x, [p * 1e-3 for p in self.m0_points], color='blue', linewidth=2)
        plt.grid(True)
        plt.ylabel(r"Total Mass, tons")
        plt.xlabel(r"Counterweight Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Space Elevator Mass vs Counterweight Altitude')

        plt.subplot(3, 1, 2)
        plt.plot(x, [p * 1e-3 for p in self.mtether_points], color='cyan', linewidth=2)
        plt.grid(True)
        plt.ylabel(r"Tether Mass, tons")
        plt.xlabel(r"Counterweight Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Tether Mass vs Counterweight Altitude')

        plt.subplot(3, 1, 3)
        plt.plot(x, [p * 1e-3 for p in self.mcw_points], color='magenta', linewidth=2)
        plt.grid(True)
        plt.ylabel(r"Counterweight Mass, tons")
        plt.xlabel(r"Counterweight Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Counterweight Mass vs Counterweight Altitude')

        plt.tight_layout(pad=2, w_pad=4, h_pad=1.0)
        plt.show()
        plt.close('all')

    def mplplotTetherProfile(self, climbers=True, figsize=(12, 14)):
        """Plot tether section profile, inner forces, and stress distribution."""

        if self.tether.anchor_safe == 'heavy' or climbers:
            self.pullforce = self.tether.pull_force
        else:
            self.pullforce = self.tether.anchor_force

        hcw  = (self.cw.hcw - self.planet.R) * 1e-6
        hsyn = (self.planet.HSYN - self.planet.R) * 1e-6
        s2   = self.tether.s2 * 1e6   # mm²

        plt.figure(figsize=figsize)

        ap, pp = self._arrangeTetherProfilePoints()
        plt.subplot(3, 1, 1)
        plt.grid(True)
        plt.xlim(0.0, hcw)
        plt.plot(ap, pp, color='0.3', linewidth=2)
        plt.fill(ap, pp, color='grey', alpha=0.5)
        plt.ylabel(r"Section Area, $\mathsf{\ mm^2}$")
        plt.xlabel(r"Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Tether Profile')
        plt.plot([hsyn, hsyn], [0.0, s2], color='0.3', ls='--', linewidth=2)
        cm = (self.se.y0 + self.planet.HSYN - self.planet.R) * 1e-6
        plt.plot([cm], [0.0], 'o', color='black', markersize=10)
        plt.text(cm, s2 * 0.1, "CM")
        plt.plot([hcw], [0.0], 'rs', markersize=20)
        plt.text(hcw * 0.85, s2 * 0.2, "Counterweight", color='red')

        alt_pts, force_pts, stress_pts = self._arrangeFSPoints(climbers)
        plt.subplot(3, 1, 2)
        plt.grid(True)
        plt.ylabel(r"Inner forces, $\mathsf{\ kN}$")
        plt.xlabel(r"Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Inner Forces Distribution Along The Tether')
        plt.xlim(0.0, hcw)
        plt.plot(alt_pts, force_pts, color='green', linewidth=2)
        plt.fill(alt_pts, force_pts, 'g', alpha=0.5)

        plt.subplot(3, 1, 3)
        plt.grid(True)
        plt.ylabel(r"Stress, $\mathsf{\ GPa}$")
        plt.xlabel(r"Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Stress Distribution Along The Tether')
        plt.xlim(0.0, hcw)
        plt.plot(alt_pts, stress_pts, color='orangered', linewidth=2)
        plt.fill(alt_pts, stress_pts, color='orangered', alpha=0.5)
        ts = self.tether.ts * 1e-9
        ws = (self.tether.ts / self.tether.ksafe) * 1e-9
        plt.plot([0.0, hcw], [ws, ws], 'b--', linewidth=2)
        plt.plot([0.0, hcw], [ts, ts], 'r--', linewidth=2)

        plt.tight_layout(pad=2, w_pad=4, h_pad=1.0)
        plt.show()
        plt.close('all')

    def showData(self, width=10, height=5, fontsize=10, digits=4):
        self.bottom  = height * 0.26 / 5.0
        self.top     = 1 - height * 0.26 / 5.0
        self.digits  = digits
        self.fontsize = fontsize
        plt.figure(figsize=(width, height))
        self._showdataPlanet()
        self._showdataTether()
        self._showdataCW()
        self._showdataSE()
        plt.tight_layout()
        plt.show()
        plt.close('all')

    def _showdataPlanet(self):
        dg = self.digits
        ax = plt.subplot2grid((5, 2), (0, 0), colspan=1, rowspan=2, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        hsyn  = self.planet.HSYN - self.planet.R
        T_    = (2 * math.pi / self.planet.OMEGA) / 3600.0
        ME    = float(5.97219e24)
        cellText = [
            [r"$\mathsf{R, \ km}$",          round(self.planet.R * 1e-3, dg)],
            [r"$\mathsf{M, \ M_{\oplus}}$",  round(self.planet.M / ME, dg)],
            [r"$\mathsf{ T, \ hr}$",          round(T_, dg)],
            [r"$\mathsf{h_{syn}, \ Mm}$",     round(hsyn * 1e-6, dg)],
        ]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("Planet: %s" % self.planet.name, fontsize=self.fontsize + 2)
        tb.set_fontsize(self.fontsize)

    def _showdataTether(self):
        th = self.tether
        dg = self.digits
        ax = plt.subplot2grid((5, 2), (0, 1), rowspan=5, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        cellText = [
            ["Taper. func. on lower part ",  th.tf_low],
            ["Taper. func. on higher part ", th.tf_high],
            [r"$\mathsf{\sigma_{b}, \ MPa}$", round(th.ts * 1e-6, dg)],
            [r"$\mathsf{k_{safe}}$",           th.ksafe],
            [r"$\mathsf{\rho, \ kg/{m^3}}$",  round(th.rho, dg)],
            [r"$\mathsf{m_{tether}, \ tonns}$", round(th.m_tether * 1e-3, dg)],
            [r"$\mathsf{s_{1}, \ mm^2}$",     round(th.s1 * 1e6, dg)],
            [r"$\mathsf{s_{2}, \ mm^2}$",     round(th.s2 * 1e6, dg)],
            [r"$\mathsf{s_{3}, \ mm^2}$",     round(th.s3 * 1e6, dg)],
            [r"$\mathsf{k_{syn}}$",            round(th.ksyn, dg)],
            [r"$\mathsf{k_{cw}}$",             round(th.kcw, dg)],
        ]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("Tether", fontsize=self.fontsize + 2)
        tb.set_fontsize(self.fontsize)

    def _showdataCW(self):
        dg  = self.digits
        ax  = plt.subplot2grid((5, 2), (4, 0), frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        hcw_ = self.cw.hcw - self.planet.R
        cellText = [
            [r"$\mathsf{h_{cw}, \ Mm}$",  round(hcw_ * 1e-6, dg)],
            [r"$\mathsf{m_{cw}, \ tons}$", round(self.cw.mcw * 1e-3, dg)],
        ]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("Counterweight", fontsize=self.fontsize + 2)
        tb.set_fontsize(self.fontsize)

    def _showdataSE(self):
        dg  = self.digits
        pl  = self.planet
        se  = self.se
        ax  = plt.subplot2grid((5, 2), (2, 0), colspan=1, rowspan=2, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        h0  = se.y0 + pl.HSYN - pl.R
        cellText = [
            [r"$\mathsf{n_{climbers}}$",         se.n_climbers],
            [r"$\mathsf{m_{climber}, \ tonns}$", round(se.m_climber * 1e-3, dg)],
            [r"$\mathsf{m_{0}, \ tonns}$",       round(se.m0 * 1e-3, dg)],
            [r"$\mathsf{h_{0}, \ Mm}$",          round(h0 * 1e-6, dg)],
        ]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("Space Elevator", fontsize=self.fontsize + 2)
        tb.set_fontsize(self.fontsize)
