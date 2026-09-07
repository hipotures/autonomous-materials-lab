# Compute Strategy

## First water milestone

W1-R is a small CPU problem. Begin with the thermal/property evaluator and lightweight provenance; neither dual-GPU acceptance nor DFT builds block this work. The GPU/DFT tiers below describe later atomistic workloads. Benchmark CFD separately; support and performance depend on the exact solver/build/model combination.

## 1. Initial hardware assumption

The local prototype is designed around:

- 2 × NVIDIA RTX 4090;
- multicore CPU;
- fast local SSD;
- optional university/HPC access.

The two GPUs should not be treated as a miniature replacement for a large DFT cluster. They are most valuable when used where consumer GPUs are exceptionally strong: machine learning, batched inference, structure relaxation with neural potentials, and molecular dynamics with suitable models.

## 2. Workload placement

### RTX 4090 — primary workload

Recommended local GPU work:

- crystal generative models;
- MACE/CHGNet-style surrogate inference;
- ML geometry relaxation;
- batched property prediction;
- embeddings and structural representation models;
- uncertainty ensembles;
- ML molecular dynamics;
- fine-tuning surrogate models.

### CPU — primary workload

Recommended CPU work:

- structure standardization;
- symmetry analysis;
- pymatgen/ASE preprocessing;
- database operations;
- workflow orchestration;
- lightweight DFT tests where GPU acceleration is not advantageous;
- post-processing.

### HPC escalation

Use HPC when a candidate justifies expensive validation:

- high-accuracy DFT with dense k-point meshes;
- large supercells;
- phonons with many displaced structures;
- extensive magnetic-state enumeration;
- ab-initio molecular dynamics;
- hybrid functionals;
- GW/BSE;
- large-scale finite-temperature calculations.

## 3. Why the RTX 4090 still matters

The RTX 4090 is not optimized for high-throughput FP64 in the way HPC accelerators are. That limits its attractiveness for some traditional electronic-structure kernels.

However, modern surrogate models typically run predominantly in lower precision and map very well to consumer GPUs.

For discovery, the important question is therefore not:

> How fast can the workstation run one final DFT calculation?

but:

> How many bad candidates can the workstation eliminate before a final DFT calculation is necessary?

That is where two 4090s can provide substantial value.

## 4. Suggested two-GPU operating modes

### Mode A — independent screening workers

```text
GPU 0 -> batch A ML relaxation
GPU 1 -> batch B ML relaxation
```

Best for high-throughput independent candidates.

### Mode B — generator + screener

```text
GPU 0 -> crystal generation
GPU 1 -> ML relaxation / scoring
```

Useful during continuous candidate production.

### Mode C — ensemble uncertainty

```text
GPU 0 -> surrogate model/checkpoint A
GPU 1 -> surrogate model/checkpoint B
```

Disagreement can contribute to epistemic-uncertainty estimation.

### Mode D — model training + production

```text
GPU 0 -> fine-tune next surrogate
GPU 1 -> continue inference with current surrogate
```

This supports an active-learning loop without stopping production.

## 5. Compute tiers

### Tier 0 — almost free

- composition filtering;
- duplicate detection;
- simple descriptors;
- database lookups;
- structural sanity checks.

### Tier 1 — cheap GPU

- ML energy;
- ML relaxation;
- ML forces/stress;
- uncertainty estimate.

Expected use: thousands to millions of candidate evaluations over time.

### Tier 2 — local reference

- modest Quantum ESPRESSO calculations;
- coarse DFT relaxation;
- selected static energies.

Expected use: tens to hundreds of candidates depending on system size.

### Tier 3 — HPC reference

- converged DFT;
- dense sampling;
- large supercells;
- phonons;
- advanced magnetic configurations.

Expected use: a small subset of candidates.

### Tier 4 — expensive research validation

- AIMD;
- free energies;
- melting-point workflows;
- hybrid functionals;
- GW/BSE;
- cross-code validation.

Expected use: only the strongest candidates.

## 6. Budget accounting

Every calculation should have a predicted and actual cost.

Suggested fields:

```text
estimated_gpu_seconds
actual_gpu_seconds
estimated_cpu_core_hours
actual_cpu_core_hours
peak_memory
wall_time
queue_time
backend
hardware_class
```

The controller can then optimize scientific return per unit of compute.

## 7. DFT presets

The orchestrator should expose named numerical presets rather than letting the AI invent arbitrary parameters.

Example:

```text
qe_relax_coarse
qe_relax_standard
qe_static_standard
qe_static_tight
qe_magnetic_scan
qe_phonon_seed
```

Each preset defines:

- pseudopotential family;
- plane-wave cutoff;
- k-point policy;
- convergence threshold;
- smearing;
- mixing;
- maximum SCF iterations;
- spin settings;
- allowed retry policy.

This is safer and more reproducible than free-form input generation.

## 8. Data locality

Large intermediate files should not be moved unnecessarily.

Suggested pattern:

```text
metadata -> relational database
structures -> database or compressed files
raw calculation directories -> local object/file store
summaries -> database
large HPC artifacts -> remote archive with indexed metadata
```

Only parsed results and selected artifacts need to reach the controller.

## 9. Early performance goal

The first benchmark should measure:

```text
candidates generated / hour
ML relaxations / GPU-hour
DFT validations / day
fraction rejected before DFT
best objective value per DFT job
```

A successful architecture should improve the final metric, not merely maximize GPU utilization.
