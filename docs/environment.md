# Development and Compute Environment

## Scope and validation status

The first W1-R benchmark needs only a CPU Python environment, NumPy/SciPy, typed schemas, a validated water-property backend and SQLite/files. The GPU and atomistic stack below is optional until the relevant later milestone. Do not reinstall a host or install all scientific engines merely to start the water model.

This document is a dependency survey, not a tested installation manifest. Selected upstream releases were checked in the [design review](design-review.md); full compatibility remains untested. Use milestone-specific locks and pin source tags before running build examples.

## 1. Recommended base platform

### Primary recommendation

Use **Ubuntu 24.04 LTS (x86_64)** as the initial workstation operating system.

Rationale:

- stable LTS base;
- mature NVIDIA driver packaging;
- officially supported by current CUDA 13.x;
- officially supported by NVIDIA HPC SDK 26.5;
- easy access to OpenMPI, BLAS/LAPACK, FFTW, HDF5 and build tooling;
- close enough to common university HPC Linux environments to reduce portability surprises.

Ubuntu 26.04 LTS is newer and is supported by CUDA 13.3, but NVIDIA HPC SDK 26.5 currently lists Ubuntu 22.04 and 24.04 as supported Ubuntu platforms. For this project, HPC compiler compatibility is more important than using the newest distribution.

Debian 13 is also viable and is supported by CUDA 13.3 and NVIDIA HPC SDK 26.5. Ubuntu 24.04 is selected mainly to reduce NVIDIA/HPC setup friction.

## 2. Hardware target

Initial workstation:

- 2 x NVIDIA GeForce RTX 4090;
- x86_64 multicore CPU;
- at least 128 GB system RAM recommended if the machine will also run DFT locally;
- NVMe storage strongly recommended;
- reliable cooling and power delivery for sustained multi-hour GPU workloads.

RTX 4090 is NVIDIA Ada, compute capability **8.9** (`sm_89`).

The workstation should be considered a high-throughput ML/screening node first and a small local DFT node second.

## 3. Version baseline

This table is a **proposed reference snapshot dated 2026-09-07**, not a known-good environment or a permanent lock file. Exact production versions will be pinned after installation and scientific acceptance tests; the review records which upstream releases were independently checked.

| Component | Current reference | Project role |
| --- | ---: | --- |
| Ubuntu | 24.04 LTS | host OS |
| NVIDIA driver | >= 580.65.06 for CUDA 13.x HPC SDK | GPU driver |
| CUDA Toolkit | 13.3.1 latest stable; do not force all software to use it | CUDA development |
| NVIDIA HPC SDK | 26.5 | NVFortran/OpenACC, QE GPU build |
| Python | **3.12 recommended for production**; 3.13+ reserved for compatibility testing | orchestration/ML |
| PyTorch | 2.14 | ML runtime |
| Quantum ESPRESSO | 7.6 | primary open-source DFT backend |
| CP2K | 2026.2 | MD / large-system DFT backend |
| ASE | 3.29.0 | atomistic interface |
| pymatgen | 2026.5.4 umbrella package | structures / phase analysis |
| pymatgen-core | fast-moving split core; track separately | core materials data structures |
| Phonopy | 4.4.0 | phonons |
| AiiDA | 2.9.2 | workflow/provenance candidate |
| MACE | 0.3.16 | ML interatomic potential |
| CHGNet | 0.4.2 | ML interatomic potential |
| MatterGen | 1.0.3 | candidate generation |

### Important versioning rule

“Latest” is not the same as “best production combination”.

The project should maintain:

- a **known-good production lock**;
- a **latest-testing lock**;
- source commit SHAs for locally modified dependencies;
- benchmark results for every toolchain update.

Do not upgrade the full stack in place without rerunning validation benchmarks.

## 4. NVIDIA stack

### Driver

Use a production NVIDIA Linux driver new enough for CUDA 13.x.

NVIDIA HPC SDK 26.5 documents a minimum CUDA 13.x driver of:

```text
580.65.06
```

A newer compatible production driver is acceptable.

Verify:

```bash
nvidia-smi
nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv
```

Expected compute capability for both RTX 4090 cards:

```text
8.9
```

### CUDA Toolkit

The latest stable CUDA Toolkit is currently 13.3.1.

Do not assume that every package should be compiled against the newest toolkit.

Recommended policy:

- system CUDA toolkit: latest validated CUDA 13.x;
- PyTorch: use its supported binary CUDA runtime;
- Quantum ESPRESSO GPU: use the CUDA runtime bundled/validated with NVIDIA HPC SDK;
- CP2K: pin the CUDA version used by its build environment.

This avoids unnecessary ABI/toolchain coupling.

### NVIDIA HPC SDK

Install NVIDIA HPC SDK 26.5.

It provides:

- `nvfortran`;
- `nvc`;
- `nvc++`;
- OpenACC;
- CUDA libraries/toolkits;
- HPC-X / OpenMPI stack;
- profiling and HPC development tools.

For Quantum ESPRESSO GPU builds, this is part of the required toolchain.

## 5. Base OS packages

A practical starting set on Ubuntu 24.04:

```bash
sudo apt update

sudo apt install -y \
  build-essential \
  git git-lfs \
  curl wget rsync \
  cmake ninja-build pkg-config \
  gfortran \
  openmpi-bin libopenmpi-dev \
  libopenblas-dev \
  libfftw3-dev libfftw3-mpi-dev \
  libscalapack-openmpi-dev \
  libhdf5-openmpi-dev \
  libxc-dev \
  python3-dev python3-venv \
  jq \
  htop \
  tmux \
  numactl \
  hwloc \
  nvme-cli \
  smartmontools
```

Not every package is required by every code. The purpose of this layer is to provide a useful host development baseline.

For CP2K, prefer its current Toolchain/Spack-assisted dependency management over manually satisfying every optional library through APT.

## 6. Python environment strategy

Do **not** install every Python package into one environment.

Use `uv` and maintain separate environments.

Suggested layout:

```text
.venv-core/
    project orchestration and schemas
    NumPy / SciPy
    validated water property backend
    SQLite / reporting

.venv-atomistic/
    ASE / pymatgen / Phonopy
    AiiDA client/plugins when needed

.venv-ml/
    PyTorch 2.14
    MACE
    CHGNet
    project ML adapters

.venv-mattergen/
    MatterGen-specific dependency set
```

### Why MatterGen should be isolated

MatterGen 1.0.3 currently pins an older stack including:

- Python 3.10 in its documented installation flow;
- `torch==2.2.1+cu118` on Linux;
- `numpy<2.0`;
- `ase>=3.22.1` in the v1.0.3 tag (do not infer an upper bound from another revision);
- older PyTorch Geometric binary dependencies.

That is intentionally incompatible with the modern main ML environment.

Initial policy:

1. run upstream MatterGen unmodified in its own Python 3.10 environment;
2. establish a reproducible baseline;
3. only then consider a project-maintained modernization branch for newer PyTorch/CUDA.

This is a good candidate for future source-level work.

## 7. Main Python environment

Recommended production interpreter:

```text
Python 3.12
```

Python 3.13 should be maintained only as a compatibility-testing target until the complete scientific stack is validated against it.

Reason:

- broad binary-wheel coverage across scientific Python packages;
- conservative compatibility target for compiled extensions and scientific libraries;
- supported by the current core stack, including PyTorch, CHGNet and AiiDA;
- lower ecosystem risk than standardizing immediately on Python 3.13 or newer;
- avoids coupling the project to the host distribution Python.

For this project, interpreter stability is more valuable than language-version novelty.

Example:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh

uv python install 3.12
uv venv .venv-core --python 3.12
uv venv .venv-ml --python 3.12
```

The exact package locks should be committed later as `uv.lock` files.

## 8. Quantum ESPRESSO

### Current reference

```text
Quantum ESPRESSO 7.6
```

QE requires:

- Fortran 2008 compiler;
- C compiler;
- MPI for distributed execution;
- OpenMP optionally;
- BLAS/LAPACK;
- FFT implementation;
- optional ScaLAPACK;
- optional HDF5;
- optional LibXC;
- NVIDIA HPC SDK for NVIDIA GPU builds.

### CPU development build

Prefer an out-of-source build.

Example:

```bash
git clone https://gitlab.com/QEF/q-e.git
cd q-e
git checkout qe-7.6

cmake -S . -B build/cpu -G Ninja \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_Fortran_COMPILER=mpif90

cmake --build build/cpu -j
```

For development, maintain separate directories:

```text
build/cpu-release
build/cpu-debug
build/gpu-release
build/gpu-debug
```

### RTX 4090 GPU build

QE 7.6 uses the current GPU configuration model based on `QE_GPU` / `QE_GPU_ARCHS` in CMake, and NVIDIA builds require NVIDIA HPC SDK/OpenACC.

Ada RTX 4090 target:

```text
sm_89
```

A development CMake configuration should therefore be based on:

```text
NVHPC compilers
+
QE_GPU_ARCHS=sm_89
```

The exact CMake invocation must be validated against the installed HPC SDK and MPI implementation and then stored as a project build preset.

The traditional `configure` flow remains available. A known GPU configuration pattern is:

```bash
./configure \
  --with-gpu=cuda \
  --with-cuda-runtime=<NVHPC CUDA version> \
  --with-cuda-cc=89 \
  --enable-openmp
```

Do not enable CUDA-aware MPI blindly on a dual-4090 workstation. RTX 4090 has no NVLink, so communication topology must be benchmarked rather than assumed beneficial.

## 9. CP2K

### Current reference

```text
CP2K 2026.2
```

Important change in the 2026 series:

- CP2K now uses CMake as the supported build system;
- 2026.2 introduced modern helper scripts for Toolchain- and Spack-based builds.

Minimum core dependencies include:

- C99 compiler;
- Fortran 2008 compiler;
- DBCSR;
- BLAS;
- LAPACK;
- MPI + ScaLAPACK for MPI builds.

### Recommended first build

Use the CP2K-managed build infrastructure before attempting a hand-tuned dependency stack.

For a development clone:

```bash
git clone https://github.com/cp2k/cp2k.git
cd cp2k
git checkout support/v2026.2
```

Then evaluate:

```bash
./make_cp2k.sh --help
```

For GPU builds, CP2K 2026.2 supports selecting a CUDA SM code through its build helper. RTX 4090 corresponds to `sm_89`.

GPU CMake builds should use Ninja.

Keep:

```text
build/release
build/debug
install/release
install/debug
```

and run CP2K regression tests after every locally modified build.

## 10. MACE

Current reference:

```text
MACE 0.3.16
```

MACE supports source installation:

```bash
git clone https://github.com/ACEsuit/mace.git
cd mace
uv pip install -e .
```

For this project, source/editable installation is preferred once the baseline is stable because:

- performance kernels can be profiled;
- GPU paths can be tested directly;
- model adapters can be patched locally;
- changes can be benchmarked before proposing upstream contributions.

Maintain an unmodified upstream baseline for comparison.

## 11. CHGNet

Current reference:

```text
CHGNet 0.4.2
```

Current package requirements include:

- Python >= 3.10;
- ASE;
- NumPy;
- pymatgen;
- PyTorch >= 2.4.1;
- Cython.

CHGNet can live in the main ML environment unless compatibility testing says otherwise.

## 12. MatterGen

Current upstream release:

```text
MatterGen 1.0.3
```

Use a dedicated upstream-compatible environment first:

```bash
git clone https://github.com/microsoft/mattergen.git
cd mattergen
git checkout v1.0.3

uv venv .venv --python 3.10
source .venv/bin/activate
uv pip install -e .
```

Also install Git LFS before cloning/pulling checkpoints:

```bash
sudo apt install git-lfs
git lfs install
```

The repository currently contains an open dependency-update effort toward much newer PyTorch. Until that modernization is upstream and tested, the project should not force MatterGen into the main PyTorch 2.14 environment.

## 13. ASE, pymatgen, Phonopy and AiiDA

Current reference versions:

```text
ASE       3.29.0
Phonopy   4.4.0
AiiDA     2.9.2
pymatgen  2026.5.4 umbrella release
```

pymatgen has recently split core functionality into a separately released `pymatgen-core` project, so its exact dependency relationship should be locked rather than inferred from “latest”.

AiiDA 2.9 is especially interesting because it now includes a built-in ZeroMQ broker, reducing local setup complexity compared with older RabbitMQ-only deployment patterns.

## 14. Source tree layout on the workstation

Do not mix project code, third-party source, build artifacts and scientific data.

Suggested layout:

```text
/srv/aml/
  project/
    autonomous-materials-lab/

  upstream/
    quantum-espresso/
    cp2k/
    mace/
    chgnet/
    mattergen/

  builds/
    qe/
    cp2k/

  envs/
    core/
    ml/
    mattergen/

  models/
    mace/
    chgnet/
    mattergen/

  pseudopotentials/
    qe/

  datasets/
    source/
    derived/

  calculations/
    scratch/
    completed/

  cache/

  logs/

  benchmarks/
```

Large generated data must not live inside the Git repository.

## 15. Pseudopotentials are part of the environment

For reproducibility, a DFT result must identify the exact pseudopotential files used.

Prepare a versioned pseudopotential store before the first serious experiment.

Store:

- source collection;
- file checksum;
- element;
- functional;
- valence configuration;
- cutoff recommendations;
- license;
- original URL/reference;
- project-local immutable identifier.

Never refer to a pseudopotential only as “the PBE oxygen potential”.

## 16. Build modes

Every source-built scientific code should have at least two modes.

### Release

```text
optimized
assertions appropriate for production
profiling available when possible
```

### Debug

```text
debug symbols
runtime checks
bounds checking where practical
sanitizers for C/C++ components where practical
reduced optimization
```

Optimization work must compare against the Release baseline and scientific regression tests.

## 17. Container policy

Do not make Docker the foundation of the initial DFT toolchain.

Reasons:

- GPU compiler/toolkit interaction;
- MPI behavior;
- HPC portability;
- profiling;
- source-level debugging.

Recommended use:

- native host builds for QE/CP2K and performance work;
- containers for CI, auxiliary services and reproducible Python tooling;
- evaluate Apptainer/Singularity later for university HPC deployment.

## 18. Reproducibility metadata

Every calculation should eventually record:

```text
OS release
kernel
CPU model
GPU model
GPU driver
CUDA runtime
compiler and version
MPI implementation/version
scientific code version
Git commit SHA
build flags
linked libraries
Python lock hash
model checkpoint hash
pseudopotential checksums
```

The environment is part of the scientific result.

## 19. Later GPU / atomistic machine validation

Before running the relevant GPU/atomistic workloads, run their hardware/software acceptance suite. This is not a prerequisite for CPU-only W1-R implementation.

### GPU

- verify both GPUs independently;
- run CUDA sample/diagnostic;
- test sustained load;
- record power, temperature and clocks;
- run PyTorch matmul/inference benchmark on each GPU;
- test concurrent independent workloads.

### CPU / memory

- memory stress test;
- BLAS benchmark;
- NUMA topology dump;
- sustained CPU load.

### Storage

- NVMe SMART health;
- sequential and random I/O benchmark;
- verify free-space monitoring;
- decide scratch retention policy.

### Scientific stack

- QE CPU reference calculation;
- QE GPU reference calculation;
- compare energies/forces numerically;
- CP2K regression subset;
- MACE GPU inference;
- CHGNet GPU relaxation;
- MatterGen generation smoke test;
- Phonopy smoke test.

Do not begin large-scale candidate generation until these tests are reproducible.
