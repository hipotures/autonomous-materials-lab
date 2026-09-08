# V5a: Water/Ethanol NRTL screening benchmark

V5a compares discrete liquid compositions using the same entry trajectory,
forebody, wall model and outlet temperature. It establishes a working property
pipeline. Its rankings remain provisional until caloric accuracy is checked
against independent measurements.

## Backend decision

- Pure fluids and exact composition endpoints: CoolProp 8.0.0 HEOS.
- Water/Ethanol mixtures: thermo 0.6.1, NRTL, GibbsExcessLiquid and FlashVL.
- No multicomponent CoolProp HEOS, imposed phases, stability-solver switches,
  interpolation across failures, or fallback to another mixture model.

The previous HEOS mixture path produced failed runs at intermediate compositions.
Their consumed masses were partial-run diagnostics, not candidate scores. The
new provider rejects requests to use HEOS for a mixture.

`property_provider.py` owns the common interface and pure-fluid provider.
`thermo_provider.py` owns the bounded binary model. Legacy `coolprop_name` configs
continue to use the pure-fluid path without importing thermo.

```yaml
coolant:
  property_provider:
    type: thermo
    model: NRTL
    parameter_set: ddbst_p05_01b
    components: [Water, Ethanol]
    composition_basis: mole
    fractions: [0.5, 0.5]
```

Mass fractions are also supported and converted to mole fractions before flash.
Zero fractions collapse to true pure fluids; small nonzero fractions stay mixtures.
The internal NRTL parameter order is **Ethanol, Water**, independent of input order.

## Thermodynamics and what has actually been verified

The parameter set reproduces the [thermo NRTL / DDBST P05.01b example](https://thermo.readthedocs.io/thermo.nrtl.html).
It uses tau_ij = B_ij/T with the documented calorie-to-joule conversion and
alpha_12 = alpha_21 = 0.2974. At 343.15 K and ethanol mole fraction 0.252,
the regression targets are gamma = [1.9360516514, 1.1536630452] and
H_excess = 582.9648539 J/mol.

These are example-reproduction targets. They are **not independent experimental
validation of heat uptake**, and the parameter fit range has not been established.
A VLE fit alone cannot establish the accuracy of its temperature derivative or
excess enthalpy. The output therefore always records
`physical_validation_complete: false` for this stage.

The [liquid phase model](https://thermo.readthedocs.io/thermo.phases.html)
uses `equilibrium_basis="Psat"` and `caloric_basis="Psat"`, including the NRTL
excess-enthalpy contribution, with ideal-gas vapor and full vapor/liquid flash.
Bulk molar enthalpy includes phase fractions and is divided by mixture molar mass
in kg/mol to obtain J/kg. Species and phase-fraction balances are checked.

For every candidate, both storage and outlet enthalpies come from the **same
provider**:

```text
delta_h = h_provider(T_out, P_out, z) - h_provider(T_storage, P_storage, z)
mass_flow = heat_removal_rate / delta_h
```

A negative absolute enthalpy can be valid because its zero is conventional.
Cross-backend absolute enthalpies must never be subtracted or used to infer a
composition discontinuity. The runner compares heat uptake near both pure
endpoints with the corresponding CoolProp heat uptake instead.

The operational screening guards are 273.15 to 500 K and 100 Pa to 2 MPa.
These guards are **not a claim of physical accuracy throughout that range**.
The ideal-gas vapor approximation becomes less reliable as pressure rises.
Liquid density uses additive pure-component volumes. Selected pure-component
correlations and package versions are recorded in the manifest; `HEOS_FIT`
correlation names do not imply a multicomponent HEOS flash.

Mixture viscosity and conductivity are unavailable and reported as null.
Two-phase Cp is also null: a phase-weighted Cp is not the equilibrium derivative
through boiling. No transport values are invented. Porous flow, tank mass,
chemistry, decomposition and coolant/trajectory mass coupling are not added here.

## Run

From `experiments/entry-evaluator`, with Python 3.12:

```bash
python -m pip install -r requirements-v5a.txt
python -m unittest discover -s tests -v
python run_mixture_sweep.py --workers 16 --output-dir mixture-v5a-nrtl
```

Use a new output directory for each run. The old `--stability-algorithm` option
has been removed. No CoolProp development build is required.

The default study uses 11 water mole fractions, storage at 293.15 K / 101325 Pa,
an outlet ceiling of 500 K, and the nominal V2 physical-study entry with the
Brandis-Johnston heating model. Chemistry is explicitly disabled.

## Outputs and failure handling

- `manifest.json`: input configs, source/config hashes, package versions,
  per-candidate model/parameter provenance and validation limitations.
- `property_validation.csv` and `property_state_grid.csv`: storage, saturation
  and the 220 requested T/P/composition states, including explicit failures.
- `continuity_report.json`: adjacent available points only; unsupported
  intervals remain gaps, and absolute cross-backend enthalpy is not compared.
- `endpoint_model_comparison.json`: near-pure delta_h versus CoolProp.
- `entry_results.csv`: consumed mass, separate `score_coolant_kg`, reference
  ratio, termination status and failure reason.
- `summary.json`: numerical completion and provisional model ranking,
  separately from physical validation.

A failed property preflight excludes the affected candidate from entry scoring
with `property_validation_failure`; supported candidates continue. The full
sweep remains incomplete. A failed reference leaves all reference ratios null.
A failed trajectory, non-finite/nonpositive mass, or explicit failure reason
cannot receive a score. Reference ratios require a successful positive-mass
reference and the same terminal condition. Skip trajectories are not eligible.
`coolant_used_kg` remains a diagnostic even on failure; use `score_coolant_kg`
for ranking. No missing composition is interpolated.

`study_complete` is retained as an alias for numerical sweep completion, not
physical validation. Exit code 0 means all requested candidates completed with
comparable numerical scores; 1 means an incomplete sweep (including excluded property candidates); 2 means
a configuration or pure-endpoint identity rejection. Numerical completion is not a flight PASS.
The existing "beats reference" count uses a 1% improvement threshold, which is
not an uncertainty estimate.

## Verification and next decision

Tests cover the documented NRTL values, excess-enthalpy inclusion, modified
Raoult-law bubble equilibrium, full two-phase balances, J/mol-to-J/kg conversion,
component ordering, mass fractions, the failed 0.4/0.5 compositions, invalid-domain
rejection, unchanged pure endpoints, coolant energy balance and failure scoring.
Single-phase Cp is checked against finite-difference enthalpy to 1e-4 relative
(the selected liquid correlations show about 6.4e-5 discrepancy at storage).

See [V5a-NRTL-results.md](V5a-NRTL-results.md) for the measured sweep results.
The next decision is whether uncertainty in delta_h can change the candidate
selection. First compare the leading candidates with an independent caloric
reference at storage and representative outlet conditions. Keep using the
current evaluator for provisional comparisons; add another mixture model or
new physical phenomena only to resolve a specific ranking ambiguity.
