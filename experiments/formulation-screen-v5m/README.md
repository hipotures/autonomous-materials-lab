# V5m-1: Known-component water formulations

The primary experimental track is now **formulations that can be prepared from
known components**, not synthesis of new molecules. This stage runs a bounded,
offline computational screen and produces a reproducible physical-test plan.
It does not operate hardware, certify handling safety, infer a contact angle,
or claim a transpiration-cooling improvement from cold-fluid properties.

## Run now: computational screen

From the repository root, use the existing V5b environment:

```bash
git pull --ff-only
cd experiments/entry-evaluator
source .venv/bin/activate
python -m unittest discover -s tests -p 'test_v5m1.py' -v
cd ../formulation-screen-v5m
python run_formulations.py --output-dir formulation-v5m1-results
```

Dependencies remain CoolProp 8.0.0, thermo 0.6.1 and chemicals 1.5.2. The existing
V5a `thermo_provider.py` is reused without changing its model or parameters.
No FeOS search, PubChem lookup or NIST crawl is run. V5e artifacts are not inputs.
No package installation or upgrade is performed by the runner.

For an interrupted run with unchanged code/configuration/environment:

```bash
python run_formulations.py --output-dir formulation-v5m1-results --resume
```

A dependency-light planning run, without numerical property calculations:

```bash
python run_formulations.py --plan-only --output-dir formulation-v5m1-plan-results
```

A plan-only directory cannot be resumed as a numerical run. Changed code,
configuration or package versions require a new output directory. Checkpoints
are content-hashed, writes atomic, and the runner holds a Linux advisory lock.
Templates are never overwritten on resume. Exit code 1 indicates failures in
supported numerical paths; unsupported models are explicit planned omissions,
not successful property predictions. Configuration/dependency errors use code 2.

## Default design: 26 formulations, all mass-based

| Additive | Mass fractions of additive | Role/model |
| --- | --- | --- |
| None | 0 | One true pure-water HEOS reference |
| Ethanol | 0.001, 0.003, 0.01, 0.03, 0.10, 0.20 | MEA cold fit; provisional V5a NRTL hot path |
| Propylene glycol | Same six fractions | MPG cold-liquid fit only |
| Glycerol | Same six fractions | MGL cold-liquid fit only |
| Sodium chloride | 0.001, 0.003, 0.01, 0.03 | MNA dissolved-solid cold-liquid fit only |
| Sodium dodecyl sulfate (SDS) | 0.0001, 0.0003, 0.001 | Measurement-only surfactant controls |

**0.001 is 0.1 mass %, not 0.001%.** Nominal 100 g recipes are computed by mass;
no volume additivity or liquid-volume percentage is assumed. Recipes assume
pure additive. Commercial concentrates, denatured alcohol, impure salts and
pre-diluted dispersions cannot be substituted without a verified composition
and a recorded preparation calculation. Batch mass is a planning setting,
not an instruction to prepare chemicals before review.

A silica-suspension catalog entry is disabled. Enabling a suspension is rejected
in this stage, rather than pretending a particle-size/pore-size ratio guarantees
no clogging. Particle distribution, aggregation, settling, exposure, filter
retention and before/after permeability require a separate protocol. Emulsions,
ternaries, immobilized beds and surface coatings also need separate designs.

## What is calculated

### Cold liquid: 293.15, 303.15 and 313.15 K at 101325 Pa

The numerical outputs are density, heat capacity, viscosity, conductivity and
enthalpy from **composition-specific aqueous-solution correlations**. They are
not weighted averages of pure-component properties. No mixture surface tension,
contact angle, solubility limit, CMC or pore-blockage probability is invented.

CoolProp's MEA, MGL and MNA fits end at 313.15 K; MPG ends at 373.15 K. Explicit
fraction/temperature checks and the library's own checks are both retained.
INCOMP models are never extended to the 500 K outlet screen. Their operational
range is not a safety limit and does not validate a boiling model.

Cold sensible heat uptake is `h(T2,p,w)-h(T1,p,w)` at identical composition and
within the same backend. The arbitrary, composition-dependent reference zeros
must not be subtracted across formulations. This sensible window is **not the
full heat uptake of an evaporating TPS coolant**.

For the same specimen, permeability and geometry, the single-phase Darcy model
implies these diagnostic ratios:

- Equal volumetric flux: pressure-drop ratio = viscosity ratio.
- Equal total mass flux: pressure-drop ratio = viscosity ratio / density ratio.

These ratios say nothing about two-phase relative permeability, changes to the
specimen, swelling, deposition or vapor blockage. They are not a TPS ranking.

### Check the water endpoint of each fitted correlation

A small additive signal can be comparable to a fit's mismatch at zero additive.
`model-endpoint-comparison.csv` therefore compares each INCOMP fit at zero
additive against HEOS water. `cold-comparison.csv` reports both ratios against
HEOS and trends relative to each fit's own zero-additive endpoint. This is a
sensitivity diagnostic, **not a bias correction or independent validation**.
There remains only one physical water reference in the formulation registry.

### Hot enthalpy window: water and water/ethanol only

The existing fixed-composition V5a NRTL flash is evaluated at 350, 400, 450 and
500 K. Both inlet and outlet enthalpy use the same provider. Bubble/dew
calculations are retained. NRTL parameter fit range and caloric accuracy remain
unvalidated; this is a provisional model comparison, not experimental evidence.
The cold MEA enthalpy is never subtracted from an NRTL outlet enthalpy.

For `r = delta_h_mixture / delta_h_water`, the heat-only mass ratio at the same
net heat load is `1/r`. At mass parity the allowed net-heat ratio is `r`; when
`r < 1`, additional net-heat reduction of `1-r` would be needed to offset that
enthalpy penalty. The code calculates **how much improvement would be needed**,
not whether wetting or a vapor layer actually provides it. No achieved surface
heat-flux suppression is inferred.

All other hot-mixture rows say `unsupported_phase_change_model`. A low-
temperature liquid fit does not provide vaporization enthalpy, a drying path,
or time-dependent vapor/liquid composition. SDS remains measurement-only.

### Retained-additive concentration stress

A separate mass balance assumes that only water leaves while all additive
remains. For initial additive mass fraction `w` and fraction `f` of initial
water removed, the residual additive fraction is:

```text
w_remaining = w / [w + (1-w)(1-f)]
```

This is explicitly a counterfactual retention stress, including for volatile
additives. It is not an evaporation, precipitation, stability or solubility
prediction. It identifies compositions that future measurements need to cover.

## Physical testing: plan first, authorization separately

The default plan has three randomized replicate rounds, with at most four
candidate formulations between a before/after pair of water controls: 117 run
slots in total. It is a proposed full screening design, not an instruction to
execute all slots immediately. Reduce the enabled catalog/fraction grid **before**
freezing a smaller campaign; do not retrospectively omit failed replicates.

The progression is preparation/identity -> storage and concentration stability
-> cold-flow compatibility -> matched thermal coupon test. Apparatus limits are
intentionally null in `protocol-template.yaml`: the actual material, pore
geometry, heater, pressure limits and available protection systems are unknown.
Models do not authorize a physical test. Local safety review, current supplier
SDS, ventilation/ignition assessment where relevant, apparatus approval and
material compatibility are prerequisites. Do not transfer numerical 500 K
outlet requests directly into physical heater settings.

Use fresh, equivalent coupons for candidate and water controls, the same material
lot/geometry/surface specification, and a verified clean fluid path. A separate
surface-conditioning experiment is needed to distinguish a beneficial deposit
from a benefit requiring additive in the flowing fluid. The current design does
not treat a modified surface as if it were an unchanged water reference.

## Later: analyze actual measurements

Copy `measurements-template.csv` to a separate measurement file and fill it from
instrument outputs or a documented export. Populate and freeze the protocol
before tests. The runner does not retrieve measurements from nonexistent sensors.
No synthetic example measurements are included in production outputs.

```bash
python analyze_bench.py \
  --plan formulation-v5m1-results/bench-plan.json \
  --measurements measured-runs.csv \
  --protocol approved-protocol.yaml \
  --output-dir bench-v5m1-results
```

Required booleans are literal `true`/`false`; blank is not a pass. `evidence_kind`
must be `measured`. A model or synthetic tag cannot become physical evidence.
Source-data, preparation and safety-review identifiers are required; this is an
audit of user-declared measurements, not independent authentication of them.

`mass_used_g` is total formulation delivered during the exposure, including
additive, with corrections documented in the preparation/data record. Incident
energy is the imposed specimen exposure: **not absorbed heat or latent heat**.
Uncertainties are declared absolute bounds; the arithmetic does not turn them
into confidence intervals. Temperature/pressure peaks are checked with their
upper bounds, candidate mass with its upper bound, and both bracketing water
masses with their lower bounds. An additional conservative incident-energy-
normalized mass ratio prevents a small exposure mismatch masquerading as savings.
Both raw mass and specific-mass comparisons must improve; no full-load scaling
law or heat-transfer coefficient is fitted from these data.

Energy, duration, hardware identity, fresh coupons, phase stability, clean flow
path, cold-flow approval, pressure, temperature, hydraulic resistance and residue
are explicit gates. Missing/failed/aborted replicates and invalid or drifting
water controls cannot be dropped to produce a winner. All planned replicates
must be valid, with at least three independent replicate blocks. A passing label
is `bench_promising_requires_confirmation`, never flight validation. Confirm any
apparent gain with a separately planned repeat campaign.

## Outputs

```text
formulation-v5m1-results/
  formulations.csv
  cold-property-grid.csv
  cold-comparison.csv
  model-endpoint-comparison.csv
  cold-sensible-window.csv
  enthalpy-window.csv
  concentration-stress.csv
  measurement-needs.csv
  bench-plan.json
  bench-plan.csv
  measurements-template.csv
  protocol-template.yaml
  summary.json
  manifest.json
  progress.json
  checkpoints/
```

The analysis command separately emits `measurement-audit.csv`,
`paired-comparisons.csv`, `summary.json` and input hashes. It never overwrites the
raw measurements or earlier molecular experiments. Missing properties stay null;
no backend failure is replaced with pure-water data.

## Tests and scientific limits

Unit tests use explicitly synthetic fixtures. Two integration tests exercise
actual CoolProp cold-mixture states and the existing NRTL path when installed.
Numerical completion is not a physical test or model accuracy guarantee. The
current code does not model transpiration, reactive cooling, tank mass, slurry
transport, dryout, local pore-scale composition, or a trajectory. No candidate
is promoted to `rankable`, and production uncertainty models are unchanged.

## Sources and existing implementation

- [CoolProp incompressible models, mass basis, fit ranges and enthalpy convention](https://coolprop.org/fluid_properties/Incompressibles.html)
- [thermo NRTL / DDBST example](https://thermo.readthedocs.io/thermo.nrtl.html)
- [Existing V5a model and limitations](../entry-evaluator/V5a.md)
- [V5a provider source](../entry-evaluator/thermo_provider.py)

Source-library coverage is not an assurance of supplier purity, formulation
stability, safe handling, or improved cooling. Supplier-specific SDS and apparatus
limits must be recorded for the actual intended physical campaign.
