# Design and Performance Evaluation of a Half-Wave Dipole Antenna in CST

Course project for TWA, B.Tech Semester 5.
Author: Anand Kumar

## Antenna Geometry

![3D Golden Half-Wave Dipole Antenna](cst/Golden%20Half-Wave%20Dipole%20Antenna.png)

*Figure 1: 3D geometry of the golden half-wave dipole antenna designed for 2.4 GHz.*

## About

This project designs a half-wave dipole for 2.4 GHz and checks its performance in two ways:

1. A CST Microwave Studio model, built with a VBA macro (`cst/`).
2. A Python method-of-moments (MoM) program that calculates the same quantities separately (`python/`).

The things measured are resonant length, S11, VSWR, input impedance, bandwidth, radiation pattern, directivity, gain, half-power beamwidth and efficiency. Results are compared with the standard theory in `docs/theory.md`.

## Design values

| Parameter | Value |
|---|---|
| Frequency | 2.4 GHz (wavelength 124.9 mm) |
| Wire radius | 1 mm |
| Feed gap | 2 mm |
| Port impedance | 50 ohm |
| Material | PEC (copper for the loss study) |
| Length | tuned to resonance, see results |

## Folder layout

```
halfwave-dipole-cst/
  cst/dipole_2p4GHz.bas    CST macro (geometry, port, monitors, mesh, solver)
  python/dipole_sim.py     MoM simulation
  docs/theory.md           theory and method description
  results/                 output of the Python run (figures, csv, summary)
  requirements.txt
  README.md
```

## How to run

Python:

```
pip install -r requirements.txt
python3 python/dipole_sim.py              # full run, takes about 10 s
python3 python/dipole_sim.py --converge   # mesh convergence check
```

The figures, `results.csv` and `summary.txt` are written to a folder called `dipole_results` (use `--out results` to write into `results/`).

CST:

1. Open CST Studio Suite and start a Microwave Studio project from the template Antennas, Wire, Time Domain.
2. Run `cst/dipole_2p4GHz.bas` from the Macros menu, or paste it into a new macro and run it.
3. In the Optimizer, vary the parameter `L` between 55 and 62 mm with the goal S11 below -20 dB at 2.4 GHz. Use parameter sweeps on `a` (radius) and `g` (gap) for the parametric study.

The macro was written without access to CST and has not been run there yet, so a command may need a small change for your CST version.

## Results from the Python simulation

Wire radius 1 mm, gap 2 mm, PEC.

| Quantity | Result |
|---|---|
| Resonant length (X = 0) | 57.1 mm (0.457 wavelength) |
| Input impedance at 2.4 GHz | 75.5 + j0 ohm |
| S11 at 2.4 GHz | -13.8 dB |
| VSWR at 2.4 GHz | 1.51 |
| -10 dB bandwidth | about 2.23 to 2.50 GHz (11.5 %) |
| Directivity | 2.14 dBi |
| E-plane half-power beamwidth | 78.5 degrees |
| Efficiency | 100 % (PEC), 99.9 % (copper) |

Radius study, with the length retuned to resonate at 2.4 GHz each time: going from a 0.25 mm to a 1 mm radius shortens the resonant length from 58.8 mm to 57.1 mm and widens the -10 dB bandwidth from 7.3 % to 11.5 %.

Textbook values for a very thin half-wave dipole are 73 + j42.5 ohm at exactly half a wavelength, 2.15 dBi directivity and a 78 degree beamwidth.

## Accuracy notes

- For a very thin wire the Python program gives 73 + j42 ohm, matching theory. The radiated power matches the power delivered at the feed, and the reactance of a short dipole is within about 8 % of theory.
- For the 1 mm wire used here, the impedance away from resonance can be off by roughly 10 to 25 %, and it changes by a few percent with the mesh. Treat the CST result as the reference. The pattern, directivity and the trends are not affected much.
- The radius study stops at 1 mm. A 2 mm radius is wider than the 2 mm feed gap, and the thin-wire model does not apply.
- No balun or feed line is modelled. The antenna is driven directly by an ideal 50 ohm port.

## Tools

Python 3 with NumPy, SciPy and Matplotlib. CST Studio Suite (Microwave Studio).
