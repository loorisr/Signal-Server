# GPU Acceleration Plan for Signal-Server

## Goal

Add a `-gpu` flag to batch-process receiver points on the GPU, accelerating the ITWOM 3.0
propagation model (and optionally LOS). The speedup comes from the **outer loop**, not from
parallelizing a single path calculation (which has sequential data dependencies).

---

## Architecture Decision

ITWOM's `lrprop()` chains `adiff()` → `alos()` → `ascat()` with data dependencies — a single
path calc is mostly sequential. The real parallelism is **across receiver locations**: each call
to `point_to_point()` for a different receiver is fully independent.

Strategy: pre-fetch terrain profiles on CPU (GDAL stays on CPU), batch-transfer to GPU, run
a CUDA kernel that executes `point_to_point()` for all receiver points in parallel, then
transfer results back.

---

## Files to Create / Modify

### New file: `models/itwom3.0.cu`
- GPU-compatible port of the ITWOM 3.0 math from `models/itwom3.0.cc`
- All functions marked `__device__`:
  - `lrprop()`, `lrprop2()`
  - `adiff()`, `adiff2()`
  - `alos()`, `alos2()`
  - `ascat()`, `saalos()`
  - `abq_alos()`
  - `qlrpfl()`, `qlrpfl2()`
  - `qlrps()`, `qlra()`
  - `hzns()`, `hzns2()`
  - `avar()`, `curve()`, `qerfi()`, `qerf()`
  - `z1sq1()`, `z1sq2()`, `qtile()`, `d1thx()`, `d1thx2()`
  - `point_to_point()`, `point_to_point_ITM()`
- One `__global__` kernel: `itwom_batch_kernel()`
  - Input: array of `GpuPropInput` structs (one per receiver point)
  - Output: array of `float` loss values
- Host-side launcher: `gpu_batch_itwom()` (called from los.cc)

### New file: `models/gpu_common.hh`
- `GpuPropInput` struct: terrain profile + all prop parameters for one receiver point
- `GpuPropOutput` struct: signal loss result
- Declaration of `gpu_batch_itwom()`

### Modified: `models/los.cc`
- In `PlotPropPath()`: when `-gpu` is active and model is ITWOM (pm==8) or LOS (pm==2),
  collect all receiver points for the current radial into a batch, call `gpu_batch_itwom()`
  once, then write results to the signal grid.
- Add `bool use_gpu` field to `PropagationRadius` struct (in `models/los.hh`).

### Modified: `main.cc`
- Parse `-gpu` flag in the argument loop (around line 1283 where other flags are parsed).
- Pass `use_gpu` into the `PropagationRadius` structs before spawning threads.

### Modified: `CMakeLists.txt`
- Add `find_package(CUDA)` or `enable_language(CUDA)`.
- Compile `models/itwom3.0.cu` with `nvcc`.
- Link `cudart`.
- Optionally add a cmake option `-DWITH_GPU=ON/OFF` to make it optional.

---

## Data Flow (GPU path)

```
CPU: GDAL tile fetch → terrain profile array
       ↓
CPU: build GpuPropInput[] for all receiver points in radial
       ↓
GPU: cudaMemcpy H→D
       ↓
GPU: itwom_batch_kernel<<<blocks, threads>>>()
     each thread handles one receiver point
       ↓
GPU: cudaMemcpy D→H  (float[] of loss values)
       ↓
CPU: write results to signal grid (same as today)
```

---

## Struct Design

```cpp
// gpu_common.hh
struct GpuPropInput {
    double elev[512];   // terrain profile (max points per path)
    int    elev_len;
    double dist;        // path distance (km)
    double hg[2];       // antenna heights (m)
    double wn;          // wave number
    double dh;          // terrain irregularity
    double ens;         // surface refractivity
    double gme;         // earth curvature
    double zgndreal, zgndimag;
    int    mdvar, klim;
    double conf, rel;
};

struct GpuPropOutput {
    float  loss_db;
    int    error_code;
};
```

---

## Kernel Design

```cuda
__global__ void itwom_batch_kernel(
    const GpuPropInput  *inputs,
    GpuPropOutput       *outputs,
    int                  n_points)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_points) return;

    // local copies of prop structs (device stack)
    prop_type   prop;
    propv_type  propv;
    propa_type  propa;

    point_to_point(inputs[idx], &prop, &propv, &propa,
                   &outputs[idx].loss_db, &outputs[idx].error_code);
}
```

Block size: 128 or 256 threads. Grid: `ceil(n_points / block_size)`.

---

## Command-line Interface

```
signalserver ... -gpu -pm 8 -o output
```

- `-gpu` enables GPU batching; ignored if no CUDA GPU is found at runtime.
- Works with `-pm 8` (ITWOM) and `-pm 2` (LOS).
- Falls back to CPU silently if CUDA device count == 0.

---

## Build System Changes (CMakeLists.txt)

```cmake
option(WITH_GPU "Enable CUDA GPU acceleration" OFF)

if(WITH_GPU)
    enable_language(CUDA)
    set(CMAKE_CUDA_ARCHITECTURES 75 86)   # Turing + Ampere; adjust as needed
    target_sources(signalserver PRIVATE models/itwom3.0.cu)
    target_compile_definitions(signalserver PRIVATE WITH_GPU)
    target_link_libraries(signalserver PRIVATE cudart)
endif()
```

When `WITH_GPU` is OFF, the `-gpu` flag at runtime prints a warning and falls back to CPU.

---

## LOS Model GPU (bonus, simpler)

`PlotLOSPath()` in `los.cc` is a simple horizon-angle sweep — much easier to port:
- Each receiver point: iterate terrain profile, track max horizon angle, check LOS.
- No complex transcendental chains.
- Port as a second kernel `los_batch_kernel()` in same `.cu` file.

---

## Open Questions (ask user before coding)

- [ ] CUDA toolkit installed? (`nvcc --version`)
- [ ] Target GPU architecture (sm_75 Turing / sm_86 Ampere / sm_89 Ada / other)?
- [ ] Max terrain profile length (currently `ARRAYSIZE` — check common.hh)?
- [ ] Prefer unified memory (`cudaMallocManaged`) or explicit transfers?
- [ ] Should `-gpu` be silently ignored or error when CUDA unavailable?
