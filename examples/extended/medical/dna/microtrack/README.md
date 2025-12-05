\page Examplemicrotrack Example microtrack

# Geant4-DNA microtrack example

Author: M. Galocha-Oliva, A. Baratto-Roldan, M.A. Cortés-Giraldo
Date: 04 November 2025
Email: miancortes@us.es

(c) The Geant4-DNA collaboration.


## 1. Introduction

The microtrack example presents a geometry to calculate microdosimetry
quantities for ion tracks in liquid water at a given energy, ensuring secondary
electron equilibrium when needed.
At each event, a random track hit defines a spherical scoring site; weighted
sampling is then used to estimate energy imparted and related quantities.

**Reference paper of this example:** 
- Front. Phys. 9 (2021) 726787

This example is provided by the Geant4-DNA collaboration.

Geant4-DNA processes and models are documented at:
http://geant4-dna.org

Any report or published results obtained using the Geant4-DNA software shall
cite the following Geant4-DNA collaboration publications:
- Med. Phys. 51 (2024) 5873–5889
- Med. Phys. 45 (2018) e722-e739
- Phys. Med. 31 (2015) 861-874
- Med. Phys. 37 (2010) 4692-4708
- Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178


## 2. Geometry

- A box-shaped World of liquid water (`G4_WATER`) centered at the origin.
- A box-shaped sensitive detector volume (`SDbox`) is created inside the world,
centered at the origin. Its X/Y dimensions match those of the world.

The size of both volumes is determined by the parameters listed below:

Configurable geometry parameters:
- `/mygeom/material <name>`: material of world/SD (default `G4_WATER`).
- `/mygeom/maxRange <value> [um|mm|cm]`: maximum electron range.
- `/mygeom/hitSelRegZ <value> [um|mm|cm]`: half-thickness of the central zone,
within the SDbox, where energy transfer points (hits) are eligible for the
random placement of the site.
- `/mygeom/hitSelRegXY <value> [um|mm|cm]`: XY size of the central zone, within
the SDbox, where energy transfer points (hits) are eligible for the random
placement of the site.
- `/mygeom/siteRadius <value> [um|mm|cm]`: radius of the site.


## 3. Particle source

- Primary generator: `G4ParticleGun` with one particle per event.
The source Z0 position is typically set at the world edge, pointing towards the
SDbox. To achieve the desired energy at the SDbox center, the user may want to
set the energy considering the expected energy loss from Z0 to the center.
- Example: protons at 20.13 MeV from z_0 = -50.99 um, to achieve 20 MeV at z=0.

Configurable beam parameters (via UI):
- `/gun/particle <name>`: e.g., `proton`, `e-`, ions available via
`G4DNAGenericIonsManager`.
- `/gun/energy <value> [Energy]`: kinetic energy compensated with (dE/dz)*Z0.
- `/beam/position/Z0 <value> [Length]`: source position along Z coordinate.

The example macro `run.mac` sets `proton` at 20.13 MeV and Z0 = -50.99 um.


## 4. Physics

- Default physics: `G4EmDNAPhysics_option2`.
- Alternative lists can be selected at pre-init using:
  - `/physics/addPhysics dna_opt1|dna_opt2|dna_opt3|dna_opt4|dna_opt5|
  dna_opt6|dna_opt7|dna_opt8`
  - `/physics/addPhysics liv` (Livermore)
  - `/physics/addPhysics penelope`
  - `/physics/addPhysics em_standard_opt4`

Production cuts can be set:
- `/run/setCutForAGivenParticle gamma <L>`
- `/run/setCutForAGivenParticle e- <L>`
- `/run/setCutForAGivenParticle e+ <L>`
- `/run/setCutForAGivenParticle proton <L>`


## 5. Scoring

All energy depositions (Hits) are recorded in the sensitive detector. At end of
each event:
- One valid hit is randomly selected within the Hit Selection Region.
- A spherical site (radius `siteRadius`) centered around a random offset within
the sphere around this hit is defined. This procedure ensures that the
spherical site sampled from a hit at the Hit Selection Region edge remains
fully contained within `SDbox`.
- For each event, we compute:
  - Edep: total energy imparted inside the site.
  - Nsel: number of hits in the Hit Selection Region.
  - Nsite: number of hits in the site.
  - Weighted and squared-weighted quantities are also histogrammed.
  - KinE_in / KinE_out at the SD boundary for primary particles.

Output produced via Geant4 analysis (ROOT by default, multi-thread merge enabled):
- File: `microtrack.root` in the run directory.
  - Histograms:
    - H1[0]: f(E_{dep}) (single event energy imparted)
    - H1[1]: E_{dep} f(E_{dep}) (weighted single event energy imparted)
    - H1[2]: E_{dep}^2 f(E_{dep}) (squared-weighted single event energy imparted)
    - H1[3]: f(y) (lineal energy)
    - H1[4]: y f(y) (weighted lineal energy)
    - H1[5]: y^2 f(y) (squared-weighted lineal energy)
    - H1[6]: f(z) (specific energy)
    - H1[7]: z f(z) (weighted specific energy)
    - H1[8]: z^2 f(z) (squared-weighted specific energy)
    - H1[9]: N_{sel} (number of hits within the central selection region)
    - H1[10]: N_{site} (number of hits within the site)
    - H1[11]: N_{int} (number of hits within the site _and_ the central selection region)
    - H1[12]: KinE_in (kinetic energy of primaries at SDbox upstream plane)
    - H1[13]: KinE_out (kinetic energy of primaries at SBbox downstream plane)
    - H2[0]: N_{site} vs E_{dep}.


## 6. Build and run

Prerequisites:
- Geant4 built with analysis, UI and visualization if desired.

Configure and build (out-of-source recommended):

```bash
mkdir build && cd build
cmake ..
make -j
```

Run in batch mode with provided macro:

```bash
./microtrack run.mac
```

Interactive with visualization:

```bash
./microtrack
```

In the UI session, the visualization is configured by `init_vis.mac` and
`vis.mac`.

Multi-threading: Set threads before initialization, e.g. in macro:
`/run/numberOfThreads 10`.


## 7. Example macro (run.mac)

The provided `run.mac` configures a water world/SD, considering
G4EmDNAPhysics_option2, proton beam at 20.13 MeV from z = -50.99 um, sites radius
of 0.5 um, and 20 events. See the file for details and as a template.


## 8. Notes and limitations

- Site selection uses a random hit within the Hit Selection Region; events
with no eligible hits are skipped from analysis.


## 9. Acknowledgments and citations

This example is provided by the Geant4-DNA collaboration.
Please cite the publications listed above.
