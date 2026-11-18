# A-SIR with segregated spaces

Python code for the aerosol-SIR (A-SIR) epidemic model with segregated indoor
spaces (outdoors, home, workplace), accompanying:

> Alexopoulos, A.H. *Aerosol Transmission and Mitigation Efficacy: Model
> Misspecification, Identifiability, and the Limits of SIR-Based Pandemic
> Policy.* Submitted for publication (2026).

It extends the single-compartment A-SIR model derived from population balance
equations in:

> Alexopoulos, A.H. (2025). Pandemic preparedness and new model developments for
> airborne diseases transmitting via aerosols. *BioSystems* 255, 105532.

## Model

A-SIR adds infected indoor spaces W (uninfected Z = 1 − W) to the SIR model:

    dS/dt = −β S I − δ S W
    dI/dt =  β S I − γ I + δ S W
    dZ/dt = −ε Z I + ζ W

βSI is near-field transmission, δSW aerosol transmission from contaminated
indoor air, εZI contamination of spaces by infected people, and ζW removal of
contamination by ventilation, deposition and viral inactivation.

The code implements a consistent model hierarchy:

| Model   | Description |
|---------|-------------|
| SIR-1   | Classical single-compartment SIR (literature reference) |
| SIR-3   | Three spaces (outdoors, home, work), near-field transmission only; reduces to SIR-1 |
| A-SIR-3 | SIR-3 plus infected-space compartments for home and work |

For each aerosol fraction (0, 20, 50, 90 %), δ is found by root-finding so that
A-SIR-3 has the same early growth rate as SIR-1, so all scenarios present the
same observable early epidemic.

## Files

| File | Purpose |
|------|---------|
| `asir3.py` | Model equations, calibration of δ, mitigation grid (Table 1 measures). Can be run on its own. |
| `paper_figs.py` | Regenerates all figures and tables of the paper (imports `asir3.py`). |
| `detect_check.py` | ΔAIC comparison of SIR and A-SIR fits for ζ = 0.05, 1.4 and 14 day⁻¹. |
| `data/` | Reference outputs, for checking a local run against the published numbers. |

`data/mitigation_metrics_winter.csv` is the winter (poor ventilation) baseline
used in the paper (ζ = 1.4 / 2.9 day⁻¹ for home / work);
`data/mitigation_metrics_proper.csv` is the proper-ventilation sensitivity case
(ζ = 14 / 29 day⁻¹). `data/paper_tables/results.json` holds every number quoted
in the text.

## Installation

Requires Python 3 (tested with 3.13).

    pip install -r requirements.txt

## Usage

Reproduce all figures and tables of the paper (about 30 s); output goes to
`paper_out/`:

    python paper_figs.py

Run the model and mitigation grid alone (winter baseline by default):

    python asir3.py
    python asir3.py --zeta-home 14 --zeta-work 29 --tag _proper

Run the model-selection (ΔAIC) check:

    python detect_check.py

All parameters (γ, ε, ζ, κ, time budget, mitigation multipliers) are set at the
top of `asir3.py`. The code is commented for readers coming from Fortran.

## Citation

If you use this code, please cite:

- The code: Alexopoulos, A.H. (2026). A-SIR with segregated spaces (v1.0).
  Zenodo. https://doi.org/10.5281/zenodo.24339219
- The model it is based on: Alexopoulos, A.H. (2025). Pandemic preparedness
  and new model developments for airborne diseases transmitting via aerosols.
  *BioSystems* 255, 105532.

A paper describing this code and its results has been submitted for
publication. Its reference will be added here once it is published.

## Licence

MIT — see `LICENSE`.

## Contact

Aleck H. Alexopoulos, CERTH, Thessaloniki, Greece — aleck@certh.gr
