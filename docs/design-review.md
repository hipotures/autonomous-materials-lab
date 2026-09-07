# Design Review — 2026-09-07

Reviewed baseline: `29921933e6b386bd1187083c3cde6fad09b94ca9`.

## Assessment

The architecture is plausible as a staged computational research platform. It is not yet an executable application, a validated thermal-protection design or evidence that a novel fluid will outperform water. The repository contains design documents only. No machine configuration, dependency lock, solver implementation or benchmark result was available for execution in this review.

The strong decisions are keeping numerical physics outside the LLM, preserving provenance, separating numerical failures from scientific outcomes, using known-water benchmarks first, and benchmarking AI against simpler policies. Preserve these decisions.

The critical path is the system evaluator. Adding more generators or faster GPUs cannot compensate for an invalid energy balance or an uncalibrated aerodynamic interaction model.

## Findings and changes

| Priority | Finding in the reviewed baseline | Required design change |
| --- | --- | --- |
| Blocking | Water-first MVP coexists with crystal-first onboarding and DFT-only reference architecture | Make W1-R CPU-only the first milestone and gate later stacks by domain |
| Blocking | Mission values are TBD; no numerical benchmark is defined | Require a complete immutable case before numerical claims; do not invent a flight mission |
| Blocking | Effective heat sink combines enthalpy and blowing without control-volume definitions | Separate fluid energy uptake, external heat-transfer modification and nozzle energy |
| Blocking | Fixed heating histories can be mistaken for coupled trajectory validation | Distinguish R component cases from C coupled cases |
| High | Reproducibility is the principal acceptance test | Add conservation, analytic limits, convergence, validity and independent-data gates |
| High | Fluid minimum can hide hardware costs and local optimizer limitations | Report best feasible result, separate loaded/consumed/dry mass and forbid unsupported system ranking |
| High | Generic ML potential could be applied to unfamiliar fluids | Gate by checkpoint chemistry, state range, observable and independent reference tests |
| High | DFT phase stability lacks an explicit cross-dataset energy compatibility rule | Use consistent references and qualify hull completeness and temperature scope |
| High | Budget validation is not yet a crash-safe execution protocol | Add atomic reservation, request identity, reconciliation and bounded retries |
| Medium | Version list can be read as a tested software combination | Mark release existence separately from installation/build/physics validation |

See the [model contract](benchmarks/model-contract.md) for the implementation-level requirements.

## Software assessment

These are role assessments from upstream documentation, not installation certifications.

| Tool | Appropriate role | Boundary / acceptance gate |
| --- | --- | --- |
| NumPy / SciPy | CPU reduced model, integration and bounded optimization | Convergence and manufactured-solution tests; optimizer success is not feasibility |
| CoolProp | Candidate water property backend for W1-R | Validate actual state range and phase; unsupported mixtures are not automatic predictions |
| Pydantic | Typed requests and persisted configuration schemas | Cross-field physics constraints and units still need explicit validators |
| SQLite + files | Initial single-host provenance and artifact store | Single-writer/transaction design, local database filesystem, checksums and tested recovery |
| Cantera | Later thermodynamics, transport and reaction kinetics | Select/validate a mechanism; does not supply the complete porous-wall/trajectory solver |
| SU2 / NEMO | Candidate later hypersonic nonequilibrium flow backend | Demonstrate required species, wall injection, conjugate heat transfer and boundary closures on a pinned build; two-phase porous delivery may need custom coupling |
| LAMMPS | Later classical or ML molecular dynamics | Force-field validity, finite-size/time convergence and property extraction are separate work |
| CP2K | Selected electronic-structure and ab-initio MD references | Does not directly deliver bulk-fluid engineering properties without a sampling workflow |
| Quantum ESPRESSO | Crystal DFT reference and selected electronic questions | Converged, compatible pseudopotentials/settings/reference phases; not the main W1 solver |
| ASE / pymatgen | Crystal structures, adapters and phase analysis | Structural deduplication tolerances and energy compatibility must be recorded |
| Phonopy | Harmonic lattice-dynamics analysis | Supercell/q-space convergence; does not prove finite-temperature or synthesis stability |
| MatterGen | Inorganic crystal candidate generation | Not a liquid/molecular generator; isolate upstream environment |
| MACE / CHGNet | Checkpoint-specific atomistic surrogate screening | Package availability does not establish validity for liquids, reactions or new thermodynamic states |
| AiiDA | Later persistent remote workflows and provenance | Evaluate when remote jobs/restarts justify it; avoid two competing workflow state stores |
| LLM | Later budgeted proposal and hypothesis layer | Must beat fixed numerical search policies at equal total cost and verified feasibility |

## Release and dependency spot checks

Upstream release records confirm the listed PyTorch 2.14.0, QE 7.6, CP2K 2026.2, MACE 0.3.16, CHGNet 0.4.2, MatterGen 1.0.3 and AiiDA 2.9.2 references. The ASE documentation identifies 3.29.0 and the Phonopy documentation identifies 4.4.0. This does not prove that they install or run together. The full environment matrix, including remaining package versions and driver/build combinations, still needs executable checks.

MatterGen's **v1.0.3 tag** declares Linux torch 2.2.1+cu118, NumPy <2 and CPython 3.10-specific binary sources for several graph extensions. It declares ASE >=3.22.1, not the ASE <=3.25.0 upper bound previously attributed to that release. A moving upstream branch may have different requirements. Pin the tag, resolver inputs, checkpoint and resolved environment together.

NVIDIA HPC SDK 26.5 includes CUDA 13.2U1 and 12.9U1 components. A system-wide newer CUDA is not a reason to force every package onto that runtime. Keep existing GPU/toolchain isolation and benchmark independent GPU jobs before splitting tightly coupled workloads.

## Missing scientific and operational contracts

### Crystal track

A hull comparison must record all competing phases, reference structures, functional, Hubbard-U policy, pseudopotentials, spin treatment, pressure and correction scheme. Do not combine raw QE energies with Materials Project/VASP energies or apply database corrections blindly. Recompute compatible references or explicitly validate a cross-method scheme. An incomplete hull can produce false stability. Negative formation energy alone is insufficient, and a 0 K harmonic result does not establish synthesis or high-temperature stability.

Keep a held-out set separated by chemical/structural family and lineage when evaluating learned models. Record possible overlap with foundation-model training data. Committee disagreement is a ranking signal until calibrated against held-out errors; agreement can still occur outside the training domain. Report coverage and residuals per regime, and retain a random exploration baseline.

### Execution and provenance

Use separate immutable candidate revisions, simulation requests, attempts, observations and decisions. A timeout/OOM/parser failure is not a physical infeasibility. Cache identity includes inputs, solver/preset version, property/checkpoint hashes and domain adapter version. Store raw artifacts with checksums and independently version the parser.

Budget transactions reserve estimated resources before dispatch and reconcile actual use afterward. Retries consume budget and have explicit attempt/time limits. After a crash, reconcile scheduler/process identity before resubmitting. A repeated request must not create duplicate jobs or charges. Use process timeouts and disk limits outside the LLM. Do not assume exactly-once execution from a status flag.

Keep untrusted literature and solver output as evidence, never executable instructions. The controller receives allowlisted tools and immutable presets; it cannot expand budgets, alter validity limits or execute arbitrary shell code. Record model identifier, prompt/tool schema versions, supplied evidence, decisions and observed cost. Use a deterministic fallback policy when the model fails.

### Data and model rights

Record the license and retrieval date for each code, checkpoint, dataset and property source separately. A source-code license does not imply redistribution rights for its weights or training data. Cache source snapshots where permitted and preserve citations with derived observations.

## Next implementation sequence

1. Choose a documented synthetic or experimental W1-R case, clearly label it, and supply all currently missing inputs.
2. Implement CPU wall/fluid balances and property-domain checks; compare independent water data and analytic limits.
3. Persist one complete result, including residuals, uncertainties, failure category and artifact checksums; verify restart/replay.
4. Add bounded coolant optimization and verify the winning design with tighter numerics.
5. Add W2-R only after chamber/nozzle conservation passes; qualify aerodynamic benefit as unknown.
6. Select/calibrate a coupled aerodynamic model before W1-C/W2-C and system-mass comparisons.
7. Add known fluids, then uncertainty-aware search and an LLM only after objective benchmarks exist.

No implementation or hardware changes are claimed by this review.

## Primary sources checked

- [ASE documentation](https://docs.ase-lib.org/install.html), [Phonopy documentation](https://phonopy.github.io/phonopy/), and [LAMMPS documentation](https://docs.lammps.org/Manual.html).
- [IAPWS-95 water formulation](https://iapws.org/technical-guidance/release/IAPWS-95).

- [CoolProp interface](https://coolprop.org/coolprop/HighLevelAPI.html) and [mixtures](https://coolprop.org/fluid_properties/Mixtures.html).
- [Cantera project](https://github.com/Cantera/cantera).
- [SU2 thermochemical nonequilibrium models](https://su2code.github.io/docs_v7/Thermochemical-Nonequilibrium/).
- [MatterGen v1.0.3 dependencies](https://github.com/microsoft/mattergen/blob/v1.0.3/pyproject.toml).
- [MACE 0.3.16](https://github.com/ACEsuit/mace/releases/tag/v0.3.16) and [CHGNet model description](https://chgnet.lbl.gov/).
- [PyTorch 2.14.0](https://github.com/pytorch/pytorch/releases/tag/v2.14.0), [QE 7.6](https://gitlab.com/QEF/q-e/-/releases/qe-7.6), [CP2K 2026.2](https://github.com/cp2k/cp2k/releases/tag/v2026.2), [CHGNet 0.4.2](https://github.com/CederGroupHub/chgnet/releases/tag/v0.4.2), [AiiDA 2.9.2](https://github.com/aiidateam/aiida-core/releases/tag/v2.9.2).
- [NVIDIA HPC SDK 26.5 release notes](https://docs.nvidia.com/hpc-sdk/archive/26.5/hpc-sdk-release-notes/index.html).
- [pymatgen energy compatibility](https://pymatgen.org/pymatgen.entries.html).
- [NASA thrust equation](https://www.grc.nasa.gov/www/k-12/BGP/rktthsum.html) and [transpiration study](https://ntrs.nasa.gov/api/citations/20000012950/downloads/20000012950.pdf).
