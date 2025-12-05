\page ExampleIAEAphsp Example IAEAphsp

# IAEAphsp — Geant4 Advanced Example

Author: M.A. Cortes-Giraldo et al.
Date: 14 Oct. 2025
Email: miancortes@us.es

IAEAphsp is an advanced Geant4 example that demonstrates **reading** and 
**writing** IAEA phase-space (IAEAphsp) files (a binary `*.IAEAphsp` with 
human-readable `*.IAEAheader`) within a minimal, configurable application.
The IAEAphsp format is defined by the IAEA Nuclear Data Section; this example
shows how to use IAEAphsp files as input source (i.e. as a primary generator)
and how to produce IAEAphsp outputs at given scoring planes.

**REFERENCE PAPER:** If you use this code, please cite:

M.A. Cortes-Giraldo et al., Int J Radiat Biol 88(1-2): 200-208 (2012)
[(doi: 10.3109/09553002.2011.627977)](https://doi.org/10.3109/09553002.2011.627977)

More information on the IAEAphsp format can be found at
[https://www-nds.iaea.org/phsp/phsp.htmlx](https://www-nds.iaea.org/phsp/phsp.htmlx)

---

## Contents

The most specific files of this example are:

  - `phsp/` directory: Contains two IAEAphsp examples, 
    1. The "test" IAEAphsp file available from the IAEAphsp project website,
    2. "PSF_example", an illustrative phsp file storing 200 particles (40 of
    each kind) recorded during 1000 original histories.
  - `iaea_phsp/` directory: Contains the files defining the IAEAphsp format
  and routines.
  - Three testing macro files:
    - `test-reader.mac`: Case in which we only read an IAEAphsp file, no 
    IAEAphsp outputs.
    - `test-writer.mac`: Regular particle gun is used, particles at specific 
    phsp planes (z const) are recorded in IAEAphsp output files.
    - `test-rw.mac`: Performs both reading and writing operations with 
    IAEAphsp files.
  - In addition, `vis.mac` macro is used for an interactive session.

---

## Geometry

Just a **box world volume** is created. You can configure world half-sizes at 
runtime via the UI commands in `/my_geom/` directory:
```
/my_geom/worldXY  <L> <unit>   # world XY half-size
/my_geom/worldZ   <L> <unit>   # world Z  half-size
```
The world material is `G4_Galactic`. No detector objects are required for the
IAEAphsp writer; the scoring planes are
**mathematical planes at constant Z positions** managed by the writer stack.

---

## Physics

This example **requires** the definition of a **reference physics list**
chosen at runtime via UI command **before initialization** (i.e., before
issuing the command `/run/initialize`). The required command is:
```
/my_phys/setList  <name>   # e.g. QGSP_BIC_HP_EMZ, QGSP_BERT_HP
```
Further, related verbosity can be controlled with:
```
/my_phys/verbose  <0|1|2>
```
If no list is set before `/run/initialize`, the application will issue a
**fatal error**.

You can also configure production cuts in the macro, as usual, e.g.:
```
/run/setCut 0.1 mm
```

---

## Particle beam (used if no phsp reader is selected)

By default, 50 MeV are shot from the center, with momentum direction parallel
to the z-axis, pointing randomly towards positive or negative direction.
This is to illustrate the recording of the incremental history number variable
(also known as `n_stat`) within the output IAEAphsp files.

Besides regular commands at `/gun/` directory, the beam can be controlled with
the following UI commands:
```
/my_beam/kinE <E> <unit>   # mean kinetic energy
/my_beam/DE   <DE> <unit>  # energy distribution half-width, flat distribution
/my_beam/X0   <X0> <unit>  # mean x-position of the beam
/my_beam/Y0   <Y0> <unit>  # mean y-position of the beam
/my_beam/Z0   <Z0> <unit>  # mean z-position of the beam
/my_beam/DX   <DX> <unit>  # x-pos distribution half width, flat distribution
/my_beam/DY   <DY> <unit>  # y-pos distribution half width, flat distribution
/my_beam/DZ   <DZ> <unit>  # z-pos distribution half width, flat distribution
```

Its verbosity can be controlled with:
```
/my_beam/verbose  <0|1|2>
```

---

## Specific commands for the IAEAphsp classes

### Activation of G4IAEAphspReader/Writer objects

The activation of IAEAphsp classes, either for reading or writing purpose,
is done via UI commands defined at a messenger class of the
ActionInitialization class, under the directory `/action/`.
This design comes from the need of passing the information safely to worker
threads and ensure MT-safe operations with the IAEAphsp routines.
```
/action/IAEAphspReader/fileName <name>  # reads from <name>.IAEA* files

/action/IAEAphspWriter/namePrefix <name>  # writes <name>[_runID].IAEA* files
/action/IAEAphspWriter/zphsp    <z_phsp> <unit>  # defines phsp plane at z-pos
```

The **G4IAEAphspReader** class only reads particle **from ONE file**.
In contrast, **more than one** zphsp values can be set to **G4IAEAphspWriter**.

### IAEAphsp Reader — controls & transforms

The G4IAEAphspReader object can be controlled with the following UI commands.
These commands are only available at `Idle` state.
Please see documentation within the class for further information.

Commands relevant for particle recycling:
```
/IAEAphspReader/recycling  <n_rec>   # Each particle is created n_rec+1 times
/IAEAphspReader/axialSymmetryX  <true|false>
/IAEAphspReader/axialSymmetryY  <true|false>
/IAEAphspReader/axialSymmetryZ  <true|false>
```

Commands relevant for simulations run in parallel reading the same IAEAphsp:
```
/IAEAphspReader/numberOfParallelRuns <n_chunk>  # No. chunks defined in file
/IAEAphspReader/parallelRun <chunk>   # Defines the piece of phsp file to read
```

Commands to mimic rotations of a linac treatment head:
```
/IAEAphspReader/collimatorAngle         <angle> <unit>
/IAEAphspReader/collimatorRotationAxis  <u_coll> <v_coll> <w_coll>
/IAEAphspReader/gantryAngle          <angle> <unit>
/IAEAphspReader/gantryRotationAxis   <u_gantry> <v_gantry> <w_gantry>
/IAEAphspReader/isocenterPosition  <x_ic> <y_ic> <z_ic> <unit>
```

Commands for custom spatial transformations of the phsp file:
```
/IAEAphspReader/rotateX    <angle> <unit>
/IAEAphspReader/rotateY    <angle> <unit>
/IAEAphspReader/rotateZ    <angle> <unit>
/IAEAphspReader/rotationOrder   <123|231|312|132|321|213>
/IAEAphspReader/translate  <x> <y> <z> <unit>
```

Verbose command:
```
/IAEAphspReader/verbose  <0|1|2>
```

---

## Build

```bash
# Configure your Geant4 path (or source your env script)
export Geant4_DIR=/path/to/geant4/lib/cmake/Geant4

# Configure & build
mkdir build && cd build
cmake -DGeant4_DIR="$Geant4_DIR" ..
cmake --build . --parallel
```
This produces the executable **`IAEAphsp`** and copies the example macros and 
phsp files into the `build/` directory.

---

## Running

### Interactive UI
```bash
./IAEAphsp
```
Launches `an interactive session (Qt/terminal depending on your Geant4 build)

### Batch (macro)
Use the provided macros from the build directory:
```bash
./IAEAphsp test-reader.mac
./IAEAphsp test-writer.mac
./IAEAphsp test-rw.mac
./IAEAphsp vis.mac
```
This example accepts an **optional** thread-count as a second argument (MT 
builds):
```bash
./IAEAphsp test-reader.mac 4
```
IAEAphsp outputs are written in the working directory.

---

## How it works (diagram)

```
       +-------------------+
       |  Run Manager      |
       |  (MT if enabled)  |
       +-------------------+
                 |
                 v
  +----------------------------------------+
  | User Initializations                   |
  | - DetectorConstruction (via /my_geom/) |
  | - PhysicsList (via /my_phys)           |
  | - ActionInitialization                 |
  +----------------------------------------+
                 |
                 v
  +--------------------------------------------------------------------+
  | PrimaryGeneratorAction                                             |
  |  - uses G4IAEAphspReader if /action/IAEAphspReader/fileName is set |
  |  - else particle gun (via /my_beam/*)                              |
  +--------------------------------------------------------------------+
                 |
                 v
  +-----------------------------------------------------------------------+
  | SteppingAction                                                        |
  |  - sends eligible tracks to                                           |
  |    G4IAEAphspWriterStack if /action/IAEAphspWriter/ commands are set  |
  +-----------------------------------------------------------------------+
                 |
                 v
      +--------------------------+
      | IAEAphspRun (per thread) |
      |  - merges at EndOfRun    |
      |  - opens/writes via      |
      |    G4IAEAphspWriter      |
      +--------------------------+
```

At end-of-run, the master consolidates thread-local stacks and the writer
produces `*.IAEAphsp`/`*.IAEAheader` files (name prefix set by
`/action/IAEAphspWriter/namePrefix`) - one pair per defined **Z plane**
(each set by `/action/IAEAphspWriter/zphsp`).

---

### Note on IAEASourceIdRegistry class

This example uses a thread-safe registry (`include/IAEASourceIdRegistry.hh`)
to coordinate the `source_ID` values passed to the IAEA routines. The goal is
to keep IDs **unique and stable** across threads and runs, avoiding surprises
from the C layer’s internal allocator.

- Readers reserve one ID per worker and reuse it on new runs; the ID is
  released when the reader is destroyed (i.e. at the end of the entire job).
- Writers reserve one ID per output plane during `OpenIAEAphspOutFiles()` and
  release them in `CloseIAEAphspOutFiles()`.

This mechanism is entirely internal; **no user commands are required**.
