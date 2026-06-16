![ClarkeLib Banner](docs/images/CLARKELIB-LOGO.png)

*With respect and gratitude to Sir Arthur C. Clarke.*

**ClarkeLib** is an open library for modelling space elevator structural parameters and loads. It lets you design a tether from anchor to counterweight, populate it with climbers, and immediately see whether the structure holds — or how much it weighs.

---

## What can be modelled

- **Tether sizing** — cross-section area at the anchor, synchronous orbit, and counterweight, driven by a working-stress constraint on your chosen material
- **Taper profiles** — exponential (constant-stress), cubic, quadratic, quartic, linear, and mixed combinations for the lower and upper tether segments
- **Counterweight** — mass and altitude, solved simultaneously with the taper so the whole structure is in equilibrium
- **Climbers** — N climbers spaced evenly along the lower tether; their downward pull is included in the structural solution
- **Force and stress distribution** — tension, cross-section area, and safety factor along the full tether length, with and without climbers
- **Counterweight altitude sweep** — scan `hcw` over a range to find the lightest viable system
- **Planets** — Earth and Mars included; parameterised by mass, radius, and rotation rate

---

## Repository layout

```
ClarkeLib/
├── clarkelib.py              # Pure Python / NumPy / SciPy module  ← start here
├── clarkelib.sage            # Original SageMath module (symbolic algebra)
├── examples/
│   ├── Space Elevator Case Study.ipynb   # Jupyter notebook walkthrough
│   └── clarkelib-test.sagews             # SageMath worksheet example
├── docs/images/
└── LICENSE.txt
```

---

## Quick start — Python

Requires Python 3.8+ with **NumPy**, **SciPy**, and **Matplotlib**.

```bash
pip install numpy scipy matplotlib
```

```python
from clarkelib import SpaceElevator

material = {'rho': 1300, 'ss': 50.0, 'ksafe': 1.3}

se = SpaceElevator(
    n_climbers   = 7,
    m_climber    = 20e3,        # kg
    anchor_force = 100e3,       # N  (anchor tension)
    material     = material,
    planet       = "Earth",
    taper_func   = ('cubic', 'exp'),
    hcw          = 50e6,        # m  (counterweight altitude above surface)
    anchor_safe  = 'light',
)

se.stats()                      # print system summary

for climber in se.climbers:     # inspect individual climbers
    print(climber)

se.plots.mplplotTetherProfile(climbers=True)   # tether cross-section & stress
se.plots.mplplotMvsHcw()                        # total mass vs counterweight altitude
```

### Taper function options

| Name | Shape | Typical use |
|---|---|---|
| `'exp'` | Constant-stress exponential | Minimum mass lower or upper tether |
| `'cubic'` | S-curve polynomial | Smooth taper, easier to manufacture |
| `'quadratic'`, `'quartic'` | Polynomial variants | Trade-off between mass and stiffness |
| `'linear'`, `'sqrt'`, `'cbrt'` | Simple profiles | Rapid prototyping / comparison |

Pass a two-element tuple `taper_func=('lower', 'upper')` to mix profiles independently for each segment.

---

## Example notebook

**`examples/Space Elevator Case Study.ipynb`** walks through a complete study:

1. Classic Earth elevator with exponential taper — baseline sizing and climber loads
2. Earth elevator with cubic lower taper — compare mass, stress profile, and equilibrium accuracy
3. Counterweight altitude sweep on both configurations — find the lightest system
4. Tether profile plots with and without climbers in transit

Open it in Jupyter Lab or any compatible environment:

```bash
jupyter lab "examples/Space Elevator Case Study.ipynb"
```

---

## Legacy SageMath module

`clarkelib.sage` is the original implementation, written for [SageMath](https://www.sagemath.org/). It uses symbolic algebra — integrals and equations are solved analytically with SageMath's CAS engine rather than numerically.

**Why keep it?**

- **Educational clarity** — symbolic expressions reveal the physics directly. The taper equation `A(h) = A₀ · exp(ρ·k/(σ) · ∫ g_eff dh)` is visible as a formula, not hidden inside a numerical solver.
- **Exact results** — symbolic integration produces closed-form tether mass expressions; useful for deriving scaling laws or checking asymptotic behaviour.
- **Historical reference** — the original design choices, naming conventions, and structural assumptions are preserved and easy to trace.

To run it, upload `clarkelib.sage` and `examples/clarkelib-test.sagews` to [CoCalc](https://cocalc.com/) (formerly SageMath Cloud) and open the worksheet.

---

## Why the Python port?

`clarkelib.py` replaces symbolic algebra with `scipy.integrate.quad` and `scipy.optimize`, making the library usable by anyone with a standard Python environment — no SageMath installation required.

- Works in Jupyter, VS Code, scripts, or CI pipelines
- Faster iteration: numerical solvers converge in milliseconds
- Easier to extend with custom taper functions or new planets
- Same class interface and parameter names as the Sage version

---

## Classes

| Class | Role |
|---|---|
| `Planet` | Planet parameters (GM, R, Ω, synchronous orbit altitude) |
| `Counterweight` | Counterweight altitude and mass |
| `Climber` | Single climber — mass, position, effective downward force |
| `Tether` | Taper functions, force/stress integrals, cross-section sizing |
| `SpaceElevator` | Top-level orchestrator — builds and solves the full system |
| `ResplotSE` | Plotting: tether profile, force distribution, mass-vs-hcw sweep |

---

## License

GPL 2 — see `LICENSE.txt`.  
Copyright 2015 Vadym Pasko. NO WARRANTY IS EXPRESSED OR IMPLIED. USE AT YOUR OWN RISK. POSSIBLE RISKS EXPLAINED IN “THE FOUNTAINS OF PARADISE” BY ARTHUR C. CLARKE

