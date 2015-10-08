
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

Use in Sage environment (http://www.sagemath.org/).
For more details see README.
"""


import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import math
import warnings

_warn_msg = {
    1: "since you've specified hcw, khcw parameter is ignored"
}

_error_msg = {

    0: "'planet' parameter must be string.",

    1: "'anchor_safe' parameter must be 'heavy' or 'light'.",

    2: "taper_func parameter must be a list "
       "or tuple of tapering fuctions names.",

    3: "Integral has not been found. "
       "Try use another tapering function.",

    4: "No solution was found. " 
       "Probably material is too weak for the selected planet."
       " \n Try using better material.",

    5: "'material' parameter must be a dictionary "
       "with 'rho', 'ss' and 'ksafe' keys with numeric values.",

    6: "Initial value of the counterweight altitude must "
       "be > then altitude of the planet's synchronous orbit.",

    7: "Final value of the counterweight altitude must "
       "be > then altitude of the planet's synchronous orbit.",

    8: "hcw must be larger then altitude of the planet's synchronous orbit",

    9: "'planet' parameter must be one of the available planet names: "

}

f0 = var('f0')          # Tether linear mass density in the anchor point
rho = var('rho')        # Tether material density
sigma = var('sigma')    # Tether material tensile strength
ksafe = var('ksafe')    # Tether safety factor
R = var('R')            # Planet radius
HSYN = var('HSYN')      # Planet synchronous orbit absolute altitude (or radius)
OMEGA = var('OMEGA')    # Planet angular velocity
GM = var('GM')          # Planet gravity constant
hcw = var('hcw')        # Counterweight absolute altitude
ksyn = var('ksyn')      # Tapering coefficient on the synchronous orbit altitude
kcw = var('kcw')        # Tapering coefficient on the counterweight altitude
h = var('h')            # Variable for altitude

# Exponential tether tapering function
_taperfunc_exp = f0*exp((GM*rho*ksafe/sigma)*(
    1/R+R^2/(2*HSYN^3)-1/h-h^2/(2*HSYN^3)))

_taperfunc_low = {
    'linear': f0*((HSYN - h)/(HSYN - R) + ksyn*(h - R)/(HSYN - R)),
    'sqrt': f0*(ksyn-1)*sqrt(h-R)/sqrt(HSYN-R)+f0,
    'cbrt': f0*(ksyn-1)*(h-R)**(1/3)/(HSYN-R)**(1/3)+f0,
    'quadratic': f0*(-ksyn+1)*(h-HSYN)^2/(R-HSYN)^2+f0*ksyn,
    'cubic': f0*(-ksyn+1)*(h-HSYN)^3/(R-HSYN)^3+f0*ksyn,
    'quartic': f0*(-ksyn+1)*(h-HSYN)^4/(R-HSYN)^4+f0*ksyn,
    'exp': _taperfunc_exp
    }
# Set of tapering functions for higher part of the tether
# (above synchronous orbit altitude)
_taperfunc_high = {
    'linear': f0*(kcw*(h - HSYN)/(hcw - HSYN) + ksyn*(hcw - h)/(hcw - HSYN)),
    'quadratic': -f0*(ksyn-kcw)*(h-HSYN)^2/(hcw-HSYN)^2+f0*ksyn,
    'cubic': -f0*(ksyn-kcw)*(h-HSYN)^3/(hcw-HSYN)^3+f0*ksyn,
    'quartic': -f0*(ksyn-kcw)*(h-HSYN)^4/(hcw-HSYN)^4+f0*ksyn,
    'exp': _taperfunc_exp
    }


class Planet(object):
    """ 
    Planet of the solar system. Currently implemented Mars and Earth.

    Parameters
    ----------
    name : str
           Name pointing to one of available planets.

    TODO: add Venus and Mercury 
    """
    G = float(6.67384e-11) # Gravitational constant (m^3/kg/s^2)
    names = ("Earth", "Mars")

    def __init__(self, name):

        if type(name) is not str:
            raise ValueError(_error_msg[0])
        elif name not in self.names:
            planetlist = "'%s', "*(len(self.names)-1) + "'%s'."
            raise ValueError(_error_msg[9] + planetlist %self.names)

        self.name = name
        if name == "Earth":
            self._setEarth()
        elif name == 'Mars':
            self._setMars()
        self._setSynOrb()

    def _setEarth(self):
        self.M = float(5.97219e24)         # Mass of Earth (kg)
        self.R = float(6371e3)             # Radius of Earth (m)
        self.OMEGA = float(7.292e-5)       # Angular velocity of Earth (rad/s)

    def _setMars(self):
        # self.T = 2*math.pi/88775
        self.M = float(6.4185e23)          # Mass of Mars (kg)
        self.R = float(3396e3)             # Radius of Mars (m)
        self.OMEGA = float(7.077651712e-5) # Angular velocity of Mars (rad/s)
        # self.OMEGA = float(7.10159e-5)   # Angular velocity of Mars (rad/s)

    def _setSynOrb(self):
        self.GM = self.G*self.M
        # Altitude of the synchronous orbit (m)
        self.HSYN = (self.GM/self.OMEGA**2)**(1.0/3)


class Tether(object):
    """
    Space Elevator tether. 
    All attributes with _low prefix correspond to lower (below base planet's
    synchronous orbit altitude) tether part. All attributes with _high prefix
    correspond to higher (above base planet's synchronous orbit altitude)
    tether part.

    Parameters
    ----------

    rho : int or float
            Density of tether material, kg/m^3
    ss : int or float
            Specific strength of tether material, MYuri
    ksafe : int or float
            Safety factor
    planet : Planet object
            Points to the Planet object.
    cw : Counterweight object
            Points to Counterweight object.
    pull_force : int or float
            Sum of anchor force and weights of all climbers at 
            their altitudes.
    anchor_force : Additional force at the anchor point (on the surface
            of the planet), N.
    taper_func : tuple, optional
            Set of tether tapering functions for lower (below base planet's 
            synchronous orbit altitude) and higher (above base planet's 
            synchronous orbit altitude) tether parts. 
    climbers : list
            List of tuples representing climbers and their main parameters
            in form (climber_force, hclimber, topbound), where
                climber_force : float
                    Weight of the climber at the altitude hclimber, N.
                hclimber : int or float
                    Altitude of the climber, m.
                topbound : int or float
                    Altitude of the next climber, m.
    anchor_safe : str, optional
                Space Elevator computation parameter. Must be one of 'heavy' 
                or 'light'. 
    """

    def __init__(self, rho, ss, ksafe, planet, cw, pull_force, anchor_force, 
                 taper_func, climbers, anchor_safe):

        self.rho = float(rho)            # Density of tether (kg/m^3)
        self.ss = float(ss)              # Specific strength (MYuri)
        self.ts = self.ss*rho*1e6        # Tensile strength (Pa)
        self.ksafe = float(ksafe)        # Safety factor

        self.planet = planet
        self.cw = cw
        self.pull_force = pull_force
        # self.pull_force_light = anchor_force + climbers[0][0]
        self.anchor_force = anchor_force

        self.anchor_safe = anchor_safe
        if anchor_safe == 'heavy':
            anchor_pullforce = pull_force
        elif anchor_safe == 'light':
            anchor_pullforce = anchor_force + climbers[0][0]
        else:
            raise ValueError(_error_msg[1])

        self.f0 = self.rho*self.ksafe*anchor_pullforce/self.ts

        if type(taper_func) is list or tuple:
            self.tf_low = taper_func[0]
            self.tf_high = taper_func[1]
        else:
            raise ValueError(_error_msg[2])

        self.solveTether()

    def _setTaperFunctions(self):

        tflow_names = _taperfunc_low.keys()
        tfhigh_names = _taperfunc_high.keys()

        if type(self.tf_low) is str and type(self.tf_high) is str:
            try:
                self._taper_low = _taperfunc_low[self.tf_low]
            except:
                tfnum = len(tflow_names)-1
                tflist_ = '%s, '*tfnum  + '%s'
                tflist = tflist_ % tflow_names
                erroer_msg = ("tapering function %s "
                              "has not been found. try one of: " %self.tf_low)
                raise ValueError(erroer_msg + tflist)
            try:
                self._taper_high = _taperfunc_high[self.tf_high]
            except:
                tfnum = len(tfhigh_names)-1
                tflist_ = '%s, '*tfnum  + '%s'
                tflist = tflist_ % tfhigh_names
                erroer_msg = ("Tapering function %s "
                              "has not been found. Try one of: " %self.tf_high)
                raise ValueError(erroer_msg + tflist)
        else:
            raise ValueError(_error_msg[1])

    def _getForces(self):

        if self.tf_low != 'exp':
            self.tetherforce_low = self.force_low.subs(hcw=self.cw.hcw, 
                                                       f0=self.f0, 
                                                       ksyn=self.ksyn, 
                                                       kcw=self.kcw)
        else:
            self.tetherforce_low = self.force_low

        if self.tf_high != 'exp':
            self.tetherforce_high = self.force_high.subs(hcw=self.cw.hcw, 
                                                         f0=self.f0, 
                                                         ksyn=self.ksyn, 
                                                         kcw=self.kcw)
        else:
            self.tetherforce_high = self.force_high

    def _findForces(self, tether_parts = 'both', ksyn_fix=1):

        self.gravlaw = GM/h^2
        self.centlaw = OMEGA^2*h

        pl = self.planet
        hcw = self.cw.hcw 

        if tether_parts == 'both' or tether_parts == 'low':
            self._findForcesLow(pl.GM, pl.HSYN, pl.R, 
                                pl.OMEGA, self.cw.hcw)
            self.force_low = self.gravforce_low - self.centforce_low

        if tether_parts == 'both' or tether_parts == 'high':
            self._findForcesHigh(pl.GM, pl.HSYN, pl.R, 
                                 pl.OMEGA, self.cw.hcw, ksyn_fix)
            self.force_high = self.gravforce_high - self.centforce_high
            self.gravforce = self.gravforce_low + self.gravforce_high
            self.centforce = self.centforce_low + self.centforce_high
  
    def _checkIntegrationState(self, f1, f2):

        try:
            f1.subs(h=1)
            f2.subs(h=1)
        except:
            raise ValueError(_error_msg[3])

    def _findForcesLow(self, pl_GM, pl_HSYN, pl_R, pl_OMEGA, cw_hcw):

        if self.tf_low != 'exp':
            assume(R>0)
            # Symbolical integration of forces 
            # caused by gravity in the lower tether part
            _gravforce_hlow = integrate(self._taper_low*self.gravlaw, h)
            # Symbolical integration of forces 
            # caused by rotation in the lower tether part
            _centforce_low = integrate(self._taper_low*self.centlaw, h)

            self._checkIntegrationState(_gravforce_hlow, _centforce_low)
        if self.tf_low != 'exp':
            self._gravforce_low = _gravforce_hlow.subs(GM=pl_GM, HSYN=pl_HSYN,
                                                       R=pl_R, OMEGA=pl_OMEGA, 
                                                       hcw=cw_hcw, f0=self.f0) 

            self._centforce_low = _centforce_low.subs(GM=pl_GM, HSYN=pl_HSYN,
                                                      R=pl_R, OMEGA=pl_OMEGA, 
                                                      hcw=cw_hcw, f0=self.f0)

            self.gravforce_low = (self._gravforce_low.subs(h=pl_HSYN) 
                                 - self._gravforce_low.subs(h=pl_R))

            self.centforce_low = (self._centforce_low.subs(h=pl_HSYN) 
                                 - self._centforce_low.subs(h=pl_R))
        else:
            gravlaw = self.gravlaw.subs(GM=pl_GM)
            centlaw = self.centlaw.subs(OMEGA=pl_OMEGA)

            self.taper_low = self._taper_low.subs(GM=pl_GM, HSYN=pl_HSYN, 
                             R=pl_R, ksafe=self.ksafe, rho=self.rho, 
                             sigma=self.ts, f0=self.f0, hcw=cw_hcw)

            gf = self.taper_low*gravlaw
            cf = self.taper_low*centlaw
            self.gravforce_low = numerical_integral(gf, pl_R, pl_HSYN)[0]
            self.centforce_low = numerical_integral(cf, pl_R, pl_HSYN)[0] 

    def _findForcesHigh(self, pl_GM, pl_HSYN, pl_R, pl_OMEGA, cw_hcw, ksynfix):

        if self.tf_high != 'exp':
            assume(R>0)
            # Symbolical integration of forces 
            # caused by gravity in the lower tether part
            _gravforce_high = integrate(self._taper_high*self.gravlaw, h)
            # Symbolical integration of forces 
            # caused by rotation in the lower tether part
            _centforce_high = integrate(self._taper_high*self.centlaw, h)

            self._checkIntegrationState(_gravforce_high, _centforce_high)
        if self.tf_high != 'exp':
            self._gravforce_high = _gravforce_high.subs(GM=pl_GM, HSYN=pl_HSYN,
                                                        R=pl_R, OMEGA=pl_OMEGA, 
                                                        hcw=cw_hcw, f0=self.f0) 

            self._centforce_high = _centforce_high.subs(GM=pl_GM, HSYN=pl_HSYN, 
                                                        R=pl_R, OMEGA=pl_OMEGA, 
                                                        hcw=cw_hcw, f0=self.f0)

            self.gravforce_high = (self._gravforce_high.subs(h=cw_hcw) 
                                  - self._gravforce_high.subs(h=pl_HSYN))

            self.centforce_high = (self._centforce_high.subs(h=cw_hcw) 
                                  - self._centforce_high.subs(h=pl_HSYN))
        else:
            gravlaw = self.gravlaw.subs(GM=pl_GM)
            centlaw = self.centlaw.subs(OMEGA=pl_OMEGA)

            self.taper_high = self._taper_high.subs(GM=pl_GM, HSYN=pl_HSYN, 
                              R=pl_R, ksafe=self.ksafe, rho=self.rho, 
                              sigma=self.ts, f0=self.f0, hcw=cw_hcw)

            gf = self.taper_high*gravlaw
            cf = self.taper_high*centlaw
            self.gravforce_high = numerical_integral(gf, pl_HSYN, cw_hcw)[0]
            self.centforce_high = numerical_integral(cf, pl_HSYN, cw_hcw)[0]

    def _getTetherAreas(self, ksyn, kcw, part='both'):

        if part == 'both' or 'low':
            # Tether section area in the anchor point
            self.s1 = self.f0/self.rho
            # Tether section area in the GEO point
            self.s2 = self.f0*ksyn/self.rho
        if part == 'both' or 'high':
            # Tether section area in the CW point
            self.s3 = self.f0*kcw/self.rho

    def _findKsyn(self):

        pl = self.planet
        if self.tf_low == 'exp':
            tf_low = _taperfunc_low[self.tf_low]
            ksyn = (1/f0)*tf_low.subs(GM=pl.GM, HSYN=pl.HSYN, R=pl.R,
                                      ksafe=self.ksafe, rho=self.rho, 
                                      sigma=self.ts, h=pl.HSYN)
            self.ksyn = ksyn.n()
        else:       
            ksyn = var('ksyn')
            fsyn = self.force_low + self.pull_force
            _eqs = (fsyn)/self.s2 - self.ts/self.ksafe
            sol = solve(_eqs, ksyn, solution_dict = True)
            self.ksyn = sol[0][ksyn].n()
            if self.tf_high == 'exp':
                ksyn_exp = ((1/f0)*_taperfunc_exp.subs(GM=pl.GM, 
                                                      HSYN=pl.HSYN,
                                                      R=pl.R, 
                                                      ksafe=self.ksafe, 
                                                      rho=self.rho, 
                                                      sigma=self.ts, 
                                                      h=pl.HSYN)).n()
                self.ksyn_fix = self.ksyn/ksyn_exp

    def _findKcwMcw(self):

        tether_force = self.gravforce - self.centforce
        _eqb = tether_force - self.cw.force_cw + self.pull_force
        _eqc = self.cw.force_cw/self.s3 - self.ts/self.ksafe

        eqb = _eqb.subs(ksyn = self.ksyn)
        eqc = _eqc.subs(ksyn = self.ksyn)

        # Solve 2 equations for kcw and mcw
        sol2 = solve([eqb,eqc], mcw, kcw, solution_dict = True)  
        self.cw.mcw = sol2[0][mcw].n()
        self.kcw = sol2[0][kcw].n()

    def getTetherMass(self):

        GM = self.planet.GM
        HSYN = self.planet.HSYN
        R = self.planet.R
        hcw = self.cw.hcw
        
        taper_low = self._taper_low.subs(GM=GM, HSYN=HSYN, R=R, 
            f0=self.s1, ksyn=self.ksyn, kcw=self.kcw, hcw=hcw, 
            sigma=self.ts, ksafe=self.ksafe, rho=self.rho)

        taper_high = self._taper_high.subs(GM=GM, HSYN=HSYN, R=R, 
            f0=self.s1, ksyn=self.ksyn, kcw=self.kcw, hcw=hcw, 
            sigma=self.ts, ksafe=self.ksafe, rho=self.rho)

        volume_low = numerical_integral(taper_low, R, HSYN)[0]
        volume_high = numerical_integral(taper_high, HSYN, hcw)[0]

        self.m_tether = (volume_low + volume_high)*self.rho

    def solveTether(self):
        """ find tether basic parameters for given couterweight altitude"""

        self._setTaperFunctions()
        if self.tf_low != 'exp' and self.tf_high == 'exp':
            self._findForces(tether_parts = 'low')
            self._getTetherAreas(ksyn, kcw, part='low')
            self._findKsyn()
            self._taper_high = self._taper_high*self.ksyn_fix
            self._findForces(tether_parts = 'high')
            self._getTetherAreas(ksyn, kcw, part='high')
        else:
            self._findForces()
            self._getTetherAreas(ksyn, kcw)
            self._findKsyn()
        self._findKcwMcw()

        ksyn_bad = self.ksyn < 0
        kcw_bad = self.ksyn < 0
        f0_bad = self.f0 < 0
        mcw_bad = self.cw.mcw < 0
        if ksyn_bad or kcw_bad or f0_bad or mcw_bad:
            raise ValueError(_error_msg[4])

        self._getTetherAreas(self.ksyn, self.kcw)
        self.cw.getCounterweightForce(symbolic = False)
        self._getForces()

    def getTaperSectionLaws(self):

        ksyn = self.ksyn
        kcw = self.kcw
        s1 = self.s1
        R = self.planet.R
        HSYN = self.planet.HSYN
        hcw = self.cw.hcw
        GM = self.planet.GM
        OMEGA = self.planet.OMEGA

        self.sections_low = self._taper_low.subs(ksyn=ksyn, kcw=kcw, f0=s1, 
                            R=R, HSYN=HSYN, OMEGA=OMEGA, GM=GM, hcw=hcw, 
                            sigma=self.ts, ksafe=self.ksafe, rho=self.rho)

        self.sections_high = self._taper_high.subs(ksyn=ksyn, kcw=kcw, f0=s1, 
                             R=R, HSYN=HSYN, OMEGA=OMEGA, GM=GM, hcw=hcw, 
                             sigma=self.ts, ksafe=self.ksafe, rho=self.rho)

    def getFoceDistribution(self, pull_force):

        pl = self.planet

        if self.tf_low != 'exp':
            _gravforce_low = self._gravforce_low.subs(ksyn=self.ksyn,
                             kcw=self.kcw, f0=self.f0)

            _centforce_low = self._centforce_low.subs(ksyn=self.ksyn,
                             kcw=self.kcw, f0=self.f0)

            _force_low = _gravforce_low - _centforce_low
            force_low = _force_low.subs(h=pl.HSYN) - _force_low.subs(h=pl.R)
            # Expression for force distribution in the lower tether part
            self.fd_low = _force_low - _force_low.subs(h=pl.R) + pull_force
            # Expression for stress distribution in the lower tether part
            self.sd_low = self.fd_low/self.sections_low 

        if self.tf_high != 'exp':
            force_low = self.tetherforce_low
            _gravforce_high = self._gravforce_high.subs(ksyn=self.ksyn,
                              kcw=self.kcw, f0=self.f0)

            _centforce_high = self._centforce_high.subs(ksyn=self.ksyn,
                              kcw=self.kcw, f0=self.f0)

            _force_high = _gravforce_high - _centforce_high
            force_high = (_force_high.subs(h=self.cw.hcw) -
                          _force_high.subs(h=pl.HSYN))
            # Expression for force distribution in the higher tether part
            self.fd_high = (force_low + _force_high -
                            _force_high.subs(h=pl.HSYN) + pull_force)
            # Expression for stress distribution in the higher tether part                       
            self.sd_high = self.fd_high/self.sections_high  

    def getEquilibriumAccuracy(self):

        eqb_accuracy = (self.tetherforce_low + self.tetherforce_high + 
                        self.pull_force - self.cw.force_cw)
        return eqb_accuracy

    def getForceDistWithClimbers(self, climbers, pull_force):
        """ compose expressions for inner forces 
            along the tether lower part with climbers """

        self.climb_forcedist_low = []
        self.climb_stressdist_low = []
        self.climber_sections = []
        _forcedist_low = self.fd_low + self.anchor_force - pull_force 
        for climber in climbers:
            climber_force = climber[0]
            _forcedist_low += climber_force
            self.climb_forcedist_low.append(_forcedist_low)
            self.climb_stressdist_low.append(_forcedist_low/self.sections_low)
            self.climber_sections.append(climber[1:])


class Counterweight(object):
    """ 
    Space Elevator counterweight. 
    """
    def __init__(self, hcw, planet):

        self.planet = planet
        self.hcw = hcw
        self.getCounterweightForce()

    def getCounterweightForce(self, symbolic = True):

        pl = self.planet
        if symbolic:
            mcw = var('mcw')
        else:
            mcw = self.mcw
        # Counterweight force
        self.force_cw = (mcw*pl.OMEGA^2)*self.hcw - mcw*pl.GM/(self.hcw**2)    

        
class SpaceElevator(object):
    """ 
    Space Elevator main class.

    Evaluates all main Space Elevator design parameters for a given 
    countarweidgt altitude hcw.

    Parameters
    ----------
    n_climbers : int, optional
                Number of Space Elevator climbers that operate on the lower
                tether part (below base planet's synchronous orbit altitude).
    m_climber : int or float, optional
                Mass of a single climber, kg.
    anchor_force : Additional force at the anchor point (on the surface
                of the planet), N.
    planet : str, optional
                Name of the base planet. Should be one of the names handled
                by 'Planet' class. 
    material : dict, optional
                Main mechanical properties of the Space Elevator tether
                material.
                Structure:
                {
                    rho : int or float
                          Density of tether material, kg/m^3
                    ss : int or float
                          Specific strength of tether material, MYuri
                    ksafe : int or float
                          Safety factor
                }
    khcw : int, float or None, optional
                Ratio of counterweight altitude to the base planet's 
                synchronous orbit altitude. If khcw and hcw is None, default
                value of khcw is used, which equals 2.
    hcw : int, float or None, optional
                Altitude of counterweight, m. Altitude count starts from the 
                planet's surface.
                If hcw is None, kcw parameter is used to calculate hcw.
    taper_func : tuple, optional
                Set of tether tapering functions for lower (below base planet's 
                synchronous orbit altitude) and higher (above base planet's 
                synchronous orbit altitude) tether parts. Must be a tuple of two
                strings pointing to the names of available tapering functions,
                specified by global variables _taperfunc_low and _taperfunc_low.
    anchor_safe : str, optional
                Space Elevator computation parameter. Must be one of 'heavy' 
                or 'light'. 
                Case description:
                'heavy' - while computing Space Elevator parameters, assumed 
                          that force at the anchor point equals sum of 
                          anchor_force and all forces caused by climbers 
                          weights on their altitudes.
                'light' - while computing Space Elevator parameters, assumed 
                          that force at the anchor point equals sum of 
                          anchor_force and weight of the first climber on the
                          anchor altitude (near the planet's surface) only.
    """
    def __init__(self, n_climbers=1, m_climber=0, anchor_force=100, 
                 planet='Earth', material='default', khcw=None, hcw=None,
                 taper_func=('exp', 'exp'), anchor_safe='heavy'):

        n_climbers = self._checkNumeric('n_climbers', n_climbers, int)
        if n_climbers<1:
            self.n_climbers = 1
        else:
            self.n_climbers = n_climbers
        m_climber = self._checkNumeric('m_climber', m_climber, float)
        self.m_climber = m_climber
        anchor_force = self._checkNumeric('anchor_force', anchor_force, float)
        if khcw is not None:
            khcw = self._checkNumeric('khcw', khcw, float)
        if hcw is not None:
            hcw = self._checkNumeric('hcw', hcw, float)
        if material == 'default':
            rho = 1300
            ss = 45
            ksafe = 1.3
        else:
            try:
                rho = float(material['rho'])
                ss = float(material['ss'])
                ksafe = float(material['ksafe'])
            except:
                raise ValueError(_error_msg[5])
        
        self.planet = Planet(planet)
        self.getClimbersForce()
        pull_force = anchor_force + self.climbers_force_total

        if hcw is None:
            if khcw is not None:
                hcw = (self.planet.HSYN - self.planet.R)*khcw + self.planet.R
            else:
                hcw = (self.planet.HSYN - self.planet.R)*2 + self.planet.R
        else:
            hcw = hcw + self.planet.R
            if khcw is not None:
                warnings.warn(_warn_msg[1])
            if hcw < self.planet.HSYN:
                raise ValueError(_error_msg[8])

        # hcw = hcw + self.planet.R
        self.cw = Counterweight(hcw, self.planet)
        self.tether = Tether(rho, ss, ksafe, self.planet, self.cw, 
                             pull_force, anchor_force, taper_func, 
                             self.climbers, anchor_safe)

        self.m0, self.y0 = self._getMassCenter(self.tether, self.cw)
        self.plots = ResplotSE(self)

    def _checkNumeric(self, name, value, value_type):

        erroer_msg = 'type of %s parameter must be numeric' %name

        if value_type is float:
            try:
                v = float(value)
            except:
                raise ValueError(erroer_msg)
        elif value_type is int:
            try:
                v = int(value)
            except:
                raise ValueError(erroer_msg)
        return v

    def _getSingleClimberWeight(self, hclimber):

        climber_force = (self.m_climber*self.planet.GM/(hclimber**2) -
                        (self.m_climber*self.planet.OMEGA**2)*hclimber)

        return climber_force

    def getClimbersForce(self):

        pl = self.planet
        hspace = np.linspace(pl.R, pl.HSYN, self.n_climbers + 1)
        self.climbers_force_total = 0
        self.climbers = []
        i = 1
        for hclimber in hspace[:-1]:
            climber_force = self._getSingleClimberWeight(hclimber)
            topbound = hspace[i]
            altitude = (hclimber - pl.R)*1e-6
            self.climbers.append((climber_force, hclimber, topbound))
            self.climbers_force_total += climber_force
            i+=1

    def _getMassCenter(self, th, cw):

        th.getTetherMass()
        pl = self.planet
        m0 = th.m_tether + cw.mcw

        cw_moment = cw.mcw*(cw.hcw - pl.HSYN)
        tether_high_moment = (cw.hcw - pl.HSYN)**2*(th.s2+2*th.s3)*th.rho/6
        tether_low_moment = (pl.HSYN-pl.R)**2*(th.s2+2*th.s1)*th.rho/6
        # Total SE mass moment regarding to point on synchronous orbit
        se_moment = cw_moment + tether_high_moment - tether_low_moment
        # SE mass center offset from point on synchronous orbit
        y0 = se_moment/m0
        
        return m0, y0

    def setHcw(self, hcw, check=True):

        hcw = hcw + self.planet.R
        if check:
            hcw = self._checkNumeric('hcw', hcw, float)

        if hcw < self.planet.HSYN:
            raise ValueError(_error_msg[8])
        self.cw = Counterweight(hcw, self.planet)
        self.tether.cw = self.cw
        self.tether.solveTether()
        self.m0, self.y0 = self._getMassCenter(self.tether, self.cw)
        self.plots.update(self)

    def varyCWAltitude(self, hcw_start, hcw_end, hcw_num):
        """ find SE target parameters in a range of counterweight altitudes"""

        hcw_start = self._checkNumeric('hcw_start', hcw_start, float)
        hcw_end = self._checkNumeric('hcw_end', hcw_end, float)
        hcw_num = self._checkNumeric('hcw_num', hcw_num, int)

        if hcw_start<self.planet.HSYN:
            raise ValueError(_error_msg[6])
        if hcw_end<self.planet.HSYN:
            raise ValueError(_error_msg[7])

        self.plots.cleanPointLists()
        hspace = np.linspace(hcw_start, hcw_end, hcw_num + 1)[:-1]

        init = True
        for hcw in hspace:

            self.setHcw(hcw, check=False)
            if init:
                m0_min = self.m0
                hcw_m0min = self.cw.hcw
                init = False
            elif self.m0 < m0_min:
                m0_min = self.m0
                hcw_m0min = self.cw.hcw

            self.plots.hcw_points.append(hcw)
            self.plots.m0_points.append(self.m0)
            self.plots.mcw_points.append(self.cw.mcw)
            self.plots.mtether_points.append(self.tether.m_tether)

        self.plots.update(self)
        self.plots.hcw_max = hcw

        return hcw_m0min - self.planet.R


class ResplotSE(object):
    """ 
    Collection of plots for Space Elevator parameters analyses.

    Notes
    -----
    Altitude count on all plots starts from the planet's radius value,
    or in other words - from the planet's surface where 
    anchor point is placed.
    """    

    def __init__(self, parent):
        self.cleanPointLists()
        self.update(parent)
        
    def update(self, parent):
        self.se = parent
        self.planet = parent.planet
        self.tether = parent.tether
        self.cw = parent.cw

    def cleanPointLists(self):

        # counterweight height points list
        self.hcw_points = []
        # space elevator total mass points list
        self.m0_points = []
        # space elevator counterweight mass points list  
        self.mcw_points = []
        # space elevator tether mass points list 
        self.mtether_points = []
    
    def _getDistPoints(self, forcedist, stressdist, hstart, hend, pnum=50):
        
        altitude_points = []
        force_points = []
        stress_points = []
        hspace = np.linspace(hstart, hend, pnum)
        for h in hspace:
            force = forcedist.subs(h = h)
            stress = stressdist.subs(h = h)
            altitude_points.append(h-self.planet.R*1e-6)
            force_points.append(force)
            stress_points.append(stress)

        return altitude_points, force_points, stress_points

    def _getDistPointsExp(self, taperfunc, sections, constforce, 
                          hstart, hend, hzero, pnum=50):
        
        altitude_points = []
        force_points = []
        stress_points = []
        hspace = np.linspace(hstart, hend, pnum)
        gravlaw = self.tether.gravlaw.subs(GM=self.planet.GM)
        centlaw = self.tether.centlaw.subs(OMEGA=self.planet.OMEGA)
        for h in hspace:
            gravforce = numerical_integral(taperfunc*gravlaw, hzero, h)[0]
            centforce = numerical_integral(taperfunc*centlaw, hzero, h)[0]
            force = constforce + gravforce - centforce
            section = sections.subs(h=h)
            stress = force/section
            alt = h-self.planet.R
            altitude_points.append(alt*1e-6)
            force_points.append(force*1e-3)
            stress_points.append(stress*1e-9)

        return altitude_points, force_points, stress_points

    def _getProfilePoints(self, section_func, hstart, hend, pnum=50):
        
        altitude_points = []
        profile_points = []
        hspace = np.linspace(hstart, hend, pnum)
        for h in hspace:
            section = section_func.subs(h = h)
            altitude_points.append(h-self.planet.R*1e-6)
            profile_points.append(section)

        return altitude_points, profile_points

    def mplplotMvsHcw(self, figsize=(8,10)):
        """ plot SE components' mass on hcw dependencies """

        pl = self.planet

        hcw = (self.hcw_max - pl.R)*1e-6
        hsyn = (self.planet.HSYN - pl.R)*1e-6
        
        plt.figure(figsize=figsize)
        x = [(point - pl.R)*1e-6 for point in self.hcw_points]

        y = [point*1e-3 for point in self.m0_points]
        plt.subplot(3,1,1)
        plt.plot(x, y, color = 'blue', linewidth = 2)
        plt.grid(True)
        plt.ylabel(r"Total Mass, tons")
        plt.xlabel(r"Counterweight Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Space Elevator Mass vs Counterweight Altitude')

        y = [point*1e-3 for point in self.mtether_points]
        plt.subplot(3,1,2)
        plt.plot(x, y, color = 'cyan', linewidth = 2)
        plt.grid(True)
        plt.ylabel(r"Tether Mass, tons")
        plt.xlabel(r"Counterweight Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Tether Mass vs Counterweight Altitude')

        y = [point*1e-3 for point in self.mcw_points]
        plt.subplot(3,1,3)
        plt.plot(x, y, color = 'magenta', linewidth = 2)
        plt.grid(True)
        plt.ylabel(r"Counterweight Mass, tons")
        plt.xlabel(r"Counterweight Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Counterweight Mass vs Counterweight Altitude')

        plt.tight_layout(pad=2, w_pad=4, h_pad=1.0)
        plt.show()
        plt.close('all')

    def mplplotTetherProfile(self, climbers = True, figsize = (12,14)):
        """ plot SE tether section areas profile """

        if self.tether.anchor_safe == 'heavy' or climbers:
            self.pullforce = self.tether.pull_force
        else:
            self.pullforce = self.tether.anchor_force
        
        hcw = (self.cw.hcw - self.planet.R)*1e-6
        hsyn = (self.planet.HSYN-self.planet.R)*1e-6
        s2 = self.tether.s2*1e6
        plt.figure(figsize=figsize)

        # Plot profile
        ap, pp = self._arrangeTetherProfilePoints()

        plt.subplot(3,1,1)
        plt.grid(True)
        plt.xlim(0.0, hcw)
        plt.plot(ap, pp, color = '0.3', linewidth = 2)
        plt.fill(ap, pp, color = 'grey', alpha = 0.5)
        plt.ylabel(r"Section Area, $\mathsf{\ mm^2}$")
        plt.xlabel(r"Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Tether Profile')
        line_x = [hsyn, hsyn]
        line_y = [0., s2]
        plt.plot(line_x, line_y, color = '0.3', ls='--', linewidth = 2)
        cm = (self.se.y0 + self.planet.HSYN - self.planet.R)*1e-6
        plt.plot([cm], [0.0], 'o', color='black', markersize=10)
        plt.text(cm, s2*0.1, "CM")
        plt.plot([hcw], [0.0], 'rs', markersize=20)
        plt.text(hcw*0.85, s2*0.2, "Counterweight", color='red')

        # Plot force distribution
        alt_pts, force_pts, stress_pts = self._arrangeFSPoints(climbers)

        plt.subplot(3,1,2)
        plt.grid(True)
        plt.ylabel(r"Inner forces, $\mathsf{\ kN}$")
        plt.xlabel(r"Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Inner Forces Distribution Along The Tether')
        plt.xlim(0.0, hcw)
        plt.plot(alt_pts, force_pts, color = 'green', linewidth = 2)
        plt.fill(alt_pts, force_pts, 'g', alpha = 0.5)

        # Plot stress distribution
        plt.subplot(3,1,3)
        plt.grid(True)
        plt.ylabel(r"Stress, $\mathsf{\ GPa}$")
        plt.xlabel(r"Altitude, $\mathsf{km\times \ 10^{3}}$")
        plt.title('Stress Distribution Along The Tether')
        plt.xlim(0.0, hcw)
        plt.plot(alt_pts, stress_pts, color="orangered", linewidth = 2)
        plt.fill(alt_pts, stress_pts, color="orangered", alpha = 0.5)
        ts = self.tether.ts*1e-9
        ws = (self.tether.ts/self.tether.ksafe)*1e-9
        alt_x = [0.0, hcw]
        ws_y = [ws, ws]
        plt.plot(alt_x, ws_y, 'b--', linewidth = 2)
        ts_y = [ts, ts]
        plt.plot(alt_x, ts_y, 'r--', linewidth = 2)

        plt.tight_layout(pad=2, w_pad=4, h_pad=1.0)
        plt.show()
        plt.close('all')

    def _arrangeFSPointsLow(self, climbers):

        th = self.tether
        pl = self.planet
        se = self.se

        if th.tf_low == 'exp':
            constforce = th.anchor_force
            if climbers:
                ap_low = []
                fp_low = []
                sp_low = []
                for climber in se.climbers:
                    climbforce, hstart, hend = climber
                    constforce += climbforce
                    ap, fp, sp = self._getDistPointsExp(th.taper_low,
                                 th.sections_low, constforce, hstart, 
                                 hend, pl.R)
                    ap_low += ap
                    fp_low += fp
                    sp_low += sp
            else:
                ap_low, fp_low, sp_low = self._getDistPointsExp(th.taper_low, 
                                         th.sections_low, self.pullforce, 
                                         pl.R, pl.HSYN, pl.R)
        else:
            climbers_data = zip(th.climb_forcedist_low, 
                                th.climb_stressdist_low, 
                                th.climber_sections)
            h = var('h')
            if climbers:
                ap_low = []
                fp_low = []
                sp_low = []
                for fd, sd, sec in climbers_data:
                    fd_low = fd.subs(h = h*1e6)*1e-3
                    sd_low = sd.subs(h = h*1e6)*1e-9
                    ap, fp, sp = self._getDistPoints(fd_low, sd_low, 
                                 sec[0]*1e-6, sec[1]*1e-6)
                    ap_low += ap
                    fp_low += fp
                    sp_low += sp
            else:
                fd_low = th.fd_low.subs(h = h*1e6)*1e-3
                sd_low = th.sd_low.subs(h = h*1e6)*1e-9
                ap_low, fp_low, sp_low = self._getDistPoints(fd_low, 
                                         sd_low, pl.R*1e-6, pl.HSYN*1e-6)

        return ap_low, fp_low, sp_low

    def _arrangeFSPointsHigh(self):

        th = self.tether
        pl = self.planet
        cw = self.cw

        if th.tf_high == 'exp':
            constforce = self.pullforce + th.tetherforce_low
            ap_high, fp_high, sp_high = self._getDistPointsExp(th.taper_high, 
                                        th.sections_high, constforce, pl.HSYN, 
                                        cw.hcw, pl.HSYN)
        else:
            fd_high = th.fd_high.subs(h = h*1e6)*1e-3
            sd_high = th.sd_high.subs(h = h*1e6)*1e-9
            ap_high, fp_high, sp_high = self._getDistPoints(fd_high, sd_high, 
                                        pl.HSYN*1e-6, cw.hcw*1e-6)

        return ap_high, fp_high, sp_high

    def _arrangeFSPoints(self, climbers):

        th = self.tether
        se = self.se

        th.getFoceDistribution(self.pullforce)

        if th.tf_low != 'exp':
            th.getForceDistWithClimbers(se.climbers, self.pullforce)

        ap_low, fp_low, sp_low = self._arrangeFSPointsLow(climbers)
        ap_high, fp_high, sp_high = self._arrangeFSPointsHigh()

        altitude_points = [ap_low[0]] + ap_low + ap_high + [ap_high[-1]]
        force_points = [0.0] + fp_low + fp_high + [0.0]
        stress_points = [0.0] + sp_low + sp_high + [0.0]

        return altitude_points, force_points, stress_points

    def _arrangeTetherProfilePoints(self):

        th = self.tether
        pl = self.planet
        cw = self.cw
        
        th.getTaperSectionLaws()
        
        sections_low = th.sections_low.subs(h = h*1e6)*1e6
        sections_high = th.sections_high.subs(h = h*1e6)*1e6

        ap_low, pp_low = self._getProfilePoints(sections_low, 
                         pl.R*1e-6, pl.HSYN*1e-6)
        ap_high, pp_high = self._getProfilePoints(sections_high, 
                           pl.HSYN*1e-6, cw.hcw*1e-6)

        ap = [ap_low[0]] + ap_low + ap_high + [ap_high[-1]]
        pp = [0.0] + pp_low + pp_high + [0.0]

        return ap, pp

    def showData(self, width=10, height=5, fontsize=10, digits=4):

        
        self.bottom = height*0.26/5.0
        self.top = (1 - height*0.26/5.0)
        self.digits = digits
        self.fontsize = fontsize
        fig = plt.figure(figsize=(width, height))
        self._showdataPlanet()
        self._showdataTether()
        self._showdataCW()
        self._showdataSE()
        plt.tight_layout()
        plt.show()
        plt.close('all')

    def _showdataPlanet(self):

        dg = self.digits        
        planetname = "Planet: %s" %self.planet.name

        ax = plt.subplot2grid((5,2), (0,0), colspan=1, 
                              rowspan=2, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        hsyn = self.planet.HSYN - self.planet.R
        T_ = (2*math.pi/self.planet.OMEGA)/3600.0
        ME = float(5.97219e24)
        mplanet = self.planet.M/ME
        HSYN = [r"$\mathsf{h_{syn}, \ Mm}$",  round(hsyn*1e-6, dg)]
        R = [r"$\mathsf{R, \ km}$", round(self.planet.R*1e-3, dg)]
        M = [r"$\mathsf{M, \ M_{\oplus}}$", round(mplanet, dg)]
        T = [r"$\mathsf{ T, \ hr}$", round(T_, dg)]

        cellText = [R, M, T, HSYN]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("%s" %planetname, fontsize=self.fontsize+2)
        tb.set_fontsize(self.fontsize)

    def _showdataTether(self):

        th = self.tether
        dg = self.digits
        
        ax = plt.subplot2grid((5,2), (0, 1), rowspan=5, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        tf_low = ["Taper. func. on lower part ",  th.tf_low]
        tf_high = ["Taper. func. on higher part ",  th.tf_high]
        sigma = [r"$\mathsf{\sigma_{b}, \ MPa}$",  round(th.ts*1e-6, dg)]
        ksafe = [r"$\mathsf{k_{safe}}$",  th.ksafe]
        rho = [r"$\mathsf{\rho, \ kg/{m^3}}$",  round(th.rho, dg)]
        m_tether = [r"$\mathsf{m_{tether}, \ tonns}$",  
                    round(th.m_tether*1e-3, dg)]
        s1 = [r"$\mathsf{s_{1}, \ mm^2}$", round(th.s1*1e6, dg)]
        s2 = [r"$\mathsf{s_{2}, \ mm^2}$", round(th.s2*1e6, dg)]
        s3 = [r"$\mathsf{s_{3}, \ mm^2}$", round(th.s3*1e6, dg)]
        ksyn = [r"$\mathsf{k_{syn}}$", round(th.ksyn, dg)]
        kcw = [r"$\mathsf{k_{cw}}$", round(th.kcw, dg)]

        cellText = [tf_low, tf_high, sigma, ksafe, rho, 
                    m_tether, s1, s2, s3, ksyn, kcw]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("Tether",  fontsize=self.fontsize+2)
        tb.set_fontsize(self.fontsize)

    def _showdataCW(self):
        
        dg = self.digits 
        ax = plt.subplot2grid((5,2), (4, 0), frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        hcw_ = self.cw.hcw - self.planet.R
        hcw = [r"$\mathsf{h_{cw}, \ Mm}$",  round(hcw_*1e-6, dg)]
        mcw = [r"$\mathsf{m_{cw}, \ tons}$", round(self.cw.mcw*1e-3, dg)]

        cellText = [hcw, mcw]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("Counterweight",  fontsize=self.fontsize+2)
        tb.set_fontsize(self.fontsize)

    def _showdataSE(self):
        
        dg = self.digits
        pl = self.planet
        se = self.se
        ax = plt.subplot2grid((5,2), (2,0), colspan=1, 
                              rowspan=2, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        h0 = se.y0 + pl.HSYN - pl.R
        n_climbers = [r"$\mathsf{n_{climbers}}$",  se.n_climbers]
        m_climber = [r"$\mathsf{m_{climber}, \ tonns}$",  
                     round(se.m_climber*1e-3, dg)]
        m0 = [r"$\mathsf{m_{0}, \ tonns}$",  round(se.m0*1e-3, dg)]
        h0 = [r"$\mathsf{h_{0}, \ Mm}$", round(h0*1e-6, dg)]

        cellText = [n_climbers, m_climber, m0, h0]
        plt.subplots_adjust(bottom=self.bottom, top=self.top)
        tb = plt.table(cellText=cellText, loc='upper center', cellLoc='left')
        plt.title("Space Elevator",  fontsize=self.fontsize+2)
        tb.set_fontsize(self.fontsize)
