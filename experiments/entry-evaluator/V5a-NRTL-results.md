# V5a NRTL benchmark: executed results

This run uses the accompanying source/config SHA-256 hashes in [V5a-NRTL-results.json](V5a-NRTL-results.json).
It was executed from an edited working tree based on f901db4; the manifest records that explicitly.
The source hashes, rather than that base commit alone, identify the implementation tested.

Python 3.12; thermo 0.6.1; CoolProp 8.0.0. All 43 entry-evaluator tests passed.
All 220 property grid states were supported; all 11 candidates reached terminal velocity.
Two additional CLI integration checks verified continuation after excluding one composition
and a clean incomplete result when every candidate is excluded.
Settings: 15 km/s, -11 degrees, nose radius 1 m, cooled area 12 m2, 32 surface zones, dt 0.05 s.
Storage: 293.15 K / 101325 Pa. Outlet ceiling: 500 K. Chemistry disabled.

| Water mole fraction | Required coolant (kg) | Ratio vs water | Numerical result |
| --- | ---: | ---: | --- |
| 0.0 | 7711.1 | 2.2129 | terminal_velocity |
| 0.1 | 7304.1 | 2.0961 | terminal_velocity |
| 0.2 | 6934.6 | 1.9901 | terminal_velocity |
| 0.3 | 6555.5 | 1.8813 | terminal_velocity |
| 0.4 | 6165.2 | 1.7693 | terminal_velocity |
| 0.5 | 5762.3 | 1.6536 | terminal_velocity |
| 0.6 | 5344.8 | 1.5338 | terminal_velocity |
| 0.7 | 4910.2 | 1.4091 | terminal_velocity |
| 0.8 | 4455.7 | 1.2787 | terminal_velocity |
| 0.9 | 3977.9 | 1.1415 | terminal_velocity |
| 1.0 | 3484.6 | 1.0000 | terminal_velocity |

The model ranks pure water first on this discrete composition grid.
The previously failing x_water = 0.4 and 0.5 compositions now produce completed-run scores.
These are model predictions, not validated coolant requirements or flight qualification.

## Backend-boundary diagnostic

At the 500 K / 25 kPa outlet, near-pure thermo heat uptake differs from CoolProp by
+0.112% at the water end and +0.289% at the ethanol end. At 2 MPa the differences
grow to +3.381% and +4.661%, respectively. Exact endpoints use CoolProp.
This measures a model discrepancy, not a confidence interval or a correction factor.
The actual active-cooling pressure range in this entry is 24.897 to 293.886 kPa.
See the JSON for every comparison pressure point.

## What this allows next

Use the evaluator to compare supported candidates provisionally. Keep numerical completion
separate from caloric accuracy and system feasibility. A successful flash is insufficient
to validate excess enthalpy. The NRTL parameters reproduce a VLE documentation example;
an independent caloric reference is still needed before treating small ranking differences
as a material-selection result. No new physical phenomena or additional backend are needed
to continue using the current screening pipeline.
