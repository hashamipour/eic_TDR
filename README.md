# eic_TDR_Q2 — Clean Layout

This repository organizes your EIC DDIS **Q²** analysis into a modular, extendable structure.

## Layout

```
.
├── analysis/                # Entry points / main programs
│   ├── DDIS_Plots_Q2_OOP.cpp
│   └── DDIS_Skim_Q2.cpp
├── include/                 # Public headers exposed to others
│   ├── Plotting.hpp
│   └── Utility.hpp
├── plotting/                # Plotting implementation (links with ROOT)
│   └── Plotting.cpp
├── utility/                 # General helpers (palettes, etc.)
│   └── Utility.cpp
├── data/
│   └── filelist.txt
├── scripts/
│   └── run_ddis.sh
├── docs/
│   └── cmmnds.txt
├── CMakeLists.txt
├── Makefile                 # Wrapper around CMake
└── README.md
```

## Build Requirements

- C++17 compiler (GCC ≥ 9 or Clang ≥ 10)
- [ROOT](https://root.cern/) (with `Core`, `Hist`, `Tree`, `RIO`, `Graf`, `Gpad` components)
- CMake ≥ 3.16
- Make (GNU make)

Ensure `root-config` is on your `PATH` and `ROOT` environment is set (e.g., `source thisroot.sh`).

## Configure & Build (cmake + make)

```bash
# From repository root
make            # runs cmake -S . -B build && cmake --build build -j
```

Or manually:

```bash
cmake -S . -B build
cmake --build build -j
```

Executables produced in `build/`:

- `ddis_plots_q2`
- `ddis_skim_q2`

## Run

### Plotting executable

```bash
./build/ddis_plots_q2 /path/to/input.root
```

This reads histograms/trees from `input.root` and produces plots using ROOT.

### Skim executable

```bash
./build/ddis_skim_q2 data/filelist.txt skim_q2.root 100
```
- `data/filelist.txt` — list of input ROOT files (one per line)
- `skim_q2.root` — output skim file name
- `100` — optional max events (example placeholder)

> Adjust arguments to match your program’s `usage` string printed on invalid input.

## Include Hygiene & Modularization

- **Public headers** live in `include/`. Only put declarations needed by other translation units here.  
- **Implementation** stays in `plotting/` and `utility/` as `.cpp` files.
- `analysis/*.cpp` are the **entry points** (contain `main()`), and include only what they need:
  - `#include "Plotting.hpp"` for plotting
  - `#include "Utility.hpp"` only if using shared helpers

### Suggestions applied
- Consolidated headers into `include/` and added a small `utility` module.
- Left `SetCustomPalette()` in `plotting/Plotting.cpp` (no-arg variant) and kept extended overloads in `utility/Utility.cpp` (`SetCustomPalette(std::string)` / `(int)`). This avoids symbol conflicts and lets analysis code call either style.
- Added a static library `eicplot` to link both plotting and utility into analysis executables.
- No wildcard includes; add more headers to `include/` only if you want to expose new APIs.

## Extending

- Add new helpers in `utility/*.cpp` and declare in `include/Utility.hpp`.
- Add new plotting features in `plotting/*.cpp` with declarations in `include/*.hpp`.
- Create a new analysis program as `analysis/MyStudy.cpp` (with `main`) and add an `add_executable` in `CMakeLists.txt` mirroring the others.

## Notes

If you previously used relative includes like `#include "../utilities/Utility.hpp"`, that’s no longer necessary: the compiler now searches `include/` via CMake’s `include_directories`. Use simple includes:

```cpp
#include "Utility.hpp"
#include "Plotting.hpp"
```

TODO : **tighten include hygiene** further (remove unused headers or split large headers).
