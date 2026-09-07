# Preflight Checklist

Complete this before implementing the autonomous discovery loop.

## Host

- [ ] Ubuntu 24.04 LTS installed and fully updated
- [ ] hostname and static machine identity chosen
- [ ] time synchronization enabled
- [ ] SSH configured
- [ ] Git and Git LFS configured
- [ ] automatic unattended reboot disabled during compute jobs

## GPU

- [ ] both RTX 4090 cards visible in `nvidia-smi`
- [ ] production NVIDIA driver installed
- [ ] compute capability reported as 8.9
- [ ] sustained dual-GPU load tested
- [ ] temperatures recorded
- [ ] power draw recorded
- [ ] PCIe link width/speed verified
- [ ] GPU ordering/UUIDs recorded

## CPU and memory

- [ ] CPU topology recorded with `lscpu` / `hwloc`
- [ ] NUMA topology understood
- [ ] RAM stress test passed
- [ ] swap policy decided
- [ ] OOM behavior understood

## Storage

- [ ] NVMe SMART data healthy
- [ ] at least one dedicated calculation/scratch filesystem
- [ ] free-space alerts configured
- [ ] backup policy for source/config/results defined
- [ ] raw scratch retention policy defined
- [ ] Git repository excluded from large generated artifacts

## Build stack

- [ ] GCC/GFortran baseline works
- [ ] CMake works
- [ ] Ninja works
- [ ] OpenMPI works locally
- [ ] BLAS/LAPACK detected
- [ ] FFTW detected
- [ ] NVIDIA HPC SDK installed
- [ ] NVFortran smoke test passes

## Python

- [ ] uv installed
- [ ] Python 3.13 available
- [ ] core environment created
- [ ] ML environment created
- [ ] MatterGen legacy environment created separately
- [ ] PyTorch sees both GPUs
- [ ] each GPU can run an independent ML job

## Scientific codes

- [ ] Quantum ESPRESSO 7.6 CPU build
- [ ] Quantum ESPRESSO 7.6 GPU build
- [ ] CPU/GPU QE numerical comparison saved
- [ ] CP2K 2026.2 build
- [ ] CP2K regression subset passes
- [ ] ASE smoke test
- [ ] pymatgen smoke test
- [ ] Phonopy 4.4 smoke test
- [ ] MACE GPU smoke test
- [ ] CHGNet GPU smoke test
- [ ] MatterGen generation smoke test
- [ ] AiiDA local profile evaluated

## Scientific data

- [ ] pseudopotential collection selected
- [ ] pseudopotential checksums stored
- [ ] license/provenance recorded
- [ ] first benchmark crystal set prepared
- [ ] known reference energies/structures recorded
- [ ] units convention defined

## Reproducibility

- [ ] source commit SHA captured for every source build
- [ ] compiler version captured
- [ ] build flags captured
- [ ] linked library versions captured
- [ ] Python lock files committed
- [ ] model checkpoint hashes captured
- [ ] GPU driver/runtime captured
- [ ] benchmark outputs archived

## Operations

- [ ] long jobs survive SSH disconnect
- [ ] local job queue strategy chosen
- [ ] logs rotate safely
- [ ] disk-space guard prevents runaway jobs
- [ ] process-level timeout policy exists
- [ ] failed calculations are preserved for diagnosis
- [ ] resource accounting format defined

## Development

- [ ] benchmark corpus committed or referenced immutably
- [ ] third-party patch policy agreed
- [ ] debug/release build separation established
- [ ] profiling tools installed
- [ ] optimization results template created

## Start criterion

Do not start autonomous candidate generation until:

1. both GPUs pass sustained tests;
2. the ML stack is reproducible;
3. QE CPU and GPU runs agree within defined numerical tolerance;
4. at least one end-to-end known-material calculation can be reproduced from a clean configuration.
