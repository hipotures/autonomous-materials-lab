# Source Development, Debugging and Optimization

## 1. Principle

The project is allowed to modify upstream open-source scientific software when there is a clear reason.

Possible reasons:

- fix a reproducible bug;
- add an adapter needed by the workflow;
- improve GPU utilization;
- reduce memory use;
- remove a throughput bottleneck;
- expose structured machine-readable output;
- add instrumentation or profiling;
- improve deterministic/reproducible behavior.

Modification is not the first response to every problem. First establish an upstream baseline.

## 2. Never patch installed packages directly

Do not edit files inside:

```text
/usr/
site-packages/
random build directories
```

Instead:

```text
upstream repository
    |
    v
pinned baseline commit
    |
    v
project branch
    |
    v
tests + benchmark
    |
    v
patch / fork / upstream PR
```

Every local change must be represented by Git history.

## 3. Baseline-first workflow

Before changing a third-party code:

1. identify exact upstream release and commit;
2. reproduce the issue or performance bottleneck;
3. save a minimal input case;
4. record scientific output;
5. record timing and hardware counters;
6. create a project branch;
7. make the smallest useful change;
8. rerun correctness tests;
9. rerun performance benchmark;
10. document the result.

A speedup without a numerical comparison is not accepted.

## 4. Correctness gates

For physics codes, optimization must preserve scientific behavior within explicitly defined tolerances.

Typical comparisons:

- total energy;
- forces;
- stress;
- relaxed lattice;
- phonon frequencies;
- SCF convergence;
- iteration count;
- symmetry.

Example benchmark record:

```yaml
case: silicon-8atom-scf
baseline_commit: abc123
candidate_commit: def456

numerics:
  energy_difference_eV_atom: 2.1e-8
  max_force_difference_eV_A: 4.3e-7

performance:
  baseline_seconds: 18.42
  candidate_seconds: 13.07
  speedup: 1.409

hardware:
  gpu: RTX 4090
  driver: ...
  compiler: ...
```

## 5. Optimization priority

Optimize the pipeline in this order:

### 1. Avoid unnecessary work

Examples:

- deduplicate structures before relaxation;
- cache descriptors;
- do not rerun identical DFT inputs;
- reuse converged wavefunctions when valid;
- batch ML inference.

### 2. Improve algorithms

Examples:

- better acquisition strategy;
- better preconditioner;
- more efficient structure matching;
- fewer SCF steps;
- better batching.

### 3. Improve implementation

Examples:

- vectorization;
- memory layout;
- fewer CPU/GPU transfers;
- kernel fusion;
- asynchronous I/O;
- improved MPI decomposition.

### 4. Low-level kernel optimization

Examples:

- CUDA kernels;
- OpenACC/OpenMP offload tuning;
- cuBLAS/cuFFT usage;
- custom PyTorch extensions.

Do not start with low-level CUDA if an algorithmic change removes 90% of the work.

## 6. Profiling toolset

Prepare:

- NVIDIA Nsight Systems;
- NVIDIA Nsight Compute;
- `nvidia-smi dmon`;
- Linux `perf`;
- `time`;
- `hwloc`;
- MPI profiling where needed;
- PyTorch profiler;
- `torch.profiler`;
- CP2K timers;
- QE timing output.

For every performance issue, determine whether it is:

```text
compute-bound
memory-bandwidth-bound
PCIe-transfer-bound
CPU-bound
MPI/communication-bound
I/O-bound
scheduler/orchestration-bound
```

## 7. Dual RTX 4090 constraint

The two RTX 4090 cards have no NVLink.

Do not assume that splitting a single tightly coupled job across both GPUs will be faster.

Default strategy:

```text
GPU 0 -> independent candidate batch
GPU 1 -> independent candidate batch
```

Use multi-GPU execution for a single job only after benchmarking shows a benefit.

For ML screening, embarrassingly parallel candidate-level execution is likely to be the highest-throughput mode.

## 8. Numerical precision

The project must explicitly control precision.

Possible modes:

- FP64 for reference numerical work when required;
- FP32 for many ML inference workloads;
- mixed precision only after validation;
- no silent conversion of reference calculations to lower precision.

An optimization that changes precision is a scientific-method change, not merely an implementation optimization.

## 9. Working with code assistants

Code assistants can be used to:

- inspect unfamiliar scientific code;
- trace call graphs;
- generate tests;
- build benchmark harnesses;
- identify obvious memory copies;
- propose refactors;
- implement adapters;
- analyze compiler errors;
- inspect profiler output;
- prepare upstream-quality patches.

They should not decide that a numerical change is correct solely from code inspection.

The validation loop is:

```text
assistant proposes
      |
      v
compiler/tests
      |
      v
scientific regression
      |
      v
benchmark
      |
      v
human/research review
```

## 10. Upstream strategy

When a local change is generally useful:

1. reduce it to a clean patch;
2. add upstream-style tests;
3. follow project formatting/contribution rules;
4. open an upstream issue/PR when appropriate;
5. keep the project patch until the upstream release containing it is adopted.

Avoid long-lived private forks unless necessary.

## 11. Patch manifest

The project should maintain a machine-readable manifest such as:

```yaml
quantum-espresso:
  upstream: qe-7.6
  commit: ...
  patches:
    - id: qe-json-timing
      commit: ...
      reason: structured timing export
      upstream_pr: null

cp2k:
  upstream: v2026.2
  commit: ...
  patches: []
```

This becomes part of calculation provenance.

## 12. Benchmark corpus

Prepare a small immutable benchmark corpus before optimization work.

Suggested categories:

- tiny molecule;
- simple semiconductor;
- metal;
- magnetic crystal;
- ionic crystal;
- 20-50 atom periodic cell;
- ML relaxation batch;
- phonon test;
- I/O-heavy workflow.

Each benchmark should include expected numerical results and performance history.

## 13. CI policy

GitHub CI should handle:

- linting;
- type checking;
- unit tests;
- schema validation;
- small CPU-only scientific smoke tests.

The local workstation should handle:

- GPU tests;
- QE/CP2K builds;
- performance benchmarks.

HPC should handle:

- large regression jobs;
- scaling tests;
- expensive validation.

CI should never pretend to validate GPU/HPC behavior it did not actually run.

## 14. Optimization log

Every meaningful optimization should answer:

- what was slow;
- how it was measured;
- what changed;
- numerical difference;
- runtime difference;
- memory difference;
- whether the change scales;
- whether it is retained.

This prevents performance folklore from becoming architecture.
