\page ExampleHadr02 Example Hadr02


Example and DPMJET: 
\author V.Ivanchenko, A.Ivanchenko, \n
UrQMD: Kh Abdel-Waged et al, A. Dotti  \n
CRMC: A. Ribon (with contributions by T. Pierog and A. Tykhonov) \n
CERN, Geneva, Switzerland \n
Geant4 Associate International \n
University of Bordeaux, CENBG/IN2P3/CNRS \n
(ESA contract 22712/09/NL/AT)

This example application is providing simulation of ion beam interaction with different 
targets. Hadronic aspects of beam target interaction are demonstrated in the example 
including longitudinal profile of energy deposition, spectra of secondary  particles,
isotope production spectra. The results are presented in a form of average numbers 
and histograms. All ion/ion models of Geant4 are available. 

In addition an interface to the FORTRAN code UrQMD-1.3cr developed by Kh. Abdel-Waged et al
for the KACST/NCMP. UrQMD model by S.A.Bass et al. Prog.Part.Nucl.Phys. 41 (1998) 225
and M.Bleicher et al. J.Phys. G25 (1999) 1859.
The UrQMD physics list uses UrQMD for supported hadron-beam interactions and FTF for
hyperons. It also exposes an ion-ion path, which is not enabled by default - see the
Limitations section below.

The interface to the Cosmic Ray Monte Carlo (CRMC) allows to use generators -
such as EPOS, DPMJET, SIBYLL etc. - for hadron-nucleus and nucleus-nucleus collisions
at very high energies.

## INSTALLATION

For simulation with Geant4 native models installation procedure is the same as for 
other examples.

## HOW TO RUN

To run the example:

```
./Hadr02 <yourmacro> QGSP_BIC
```

The last parameter is optional. It is the name of Geant4 reference Physics List, 
alternatively Physics List can be defined via environment variable

```
setenv PHYSLIST QGSP_BIC
```

## ACTIVATION OF URQMD INTERFACE

The UrQMD 1.3 FORTRAN code is NOT provided with the Geant4 code-base. The interface is
written for the specific urqmd-1.3cr ("cosmic ray") version, which is request-only: it
is no longer a public download (urqmd.org now serves 3.4 / 4.0, which are not compatible
with this interface). Request urqmd-1.3cr.tar(.gz) from the UrQMD group (M. Bleicher,
ITP Frankfurt).

Unpack the tarball inside the urqmd1_3 sub-directory of this example. The archive
contains a top-level urqmd-1.3cr directory, so the sources land in urqmd1_3/urqmd-1.3cr :

```
cd examples/extended/hadronic/Hadr02/urqmd1_3
tar xf /path/to/urqmd-1.3cr.tar
```

No manual build step is needed. Enable the interface with the CMake option G4_USE_URQMD
(OFF by default); CMake then compiles the UrQMD .f sources into a static library and
links the gfortran runtime automatically. A gfortran compiler is required.

```
cd examples/extended/hadronic/Hadr02
cmake -S . -B build -DCMAKE_PREFIX_PATH=<geant4-install> -DG4_USE_URQMD=ON
cmake --build build -j
```

UrQMD can be used in two ways: only for ion-ion interactions, or as a physics list
where supported hadronic inelastic interactions use UrQMD and hyperons use FTF.

To use UrQMD only for ion-ion physics (added on top of a reference physics list).
NOTE: as shipped this ion-ion path is not usable - with the trimmed arrays a dense A+A
event overflows UrQMD's collision table and aborts (a fixable capacity limit; see the
Limitations section). Use the full UrQMD physics list (below) for hadron beams. The
ion-ion invocation is:

```
./Hadr02 urqmd.in QGSP_BIC
```

The last parameter is the reference physics list on top of which the UrQMD ion physics
is added. Alternatively the physics list can be defined via the environment variable:

```
export PHYSLIST=QGSP_BIC
```

To use the UrQMD physics list (UrQMD for supported interactions and FTF for hyperons):

```
./Hadr02 hadr02.in UrQMD
```
or:
```
export PHYSLIST=UrQMD
./Hadr02 hadr02.in
```

To reuse the UrQMD physics list in another application, copy the relevant headers and
sources (*UrQMD*) together with the urqmd1_3 sub-directory, and add the same G4_USE_URQMD
CMake block used by this example. The UrQMD interface is not thread-safe, so the reusing
application must run single-threaded (serial run manager); see the Limitations section.

### Limitations of the UrQMD interface

- Usable for hadron beams (proton, neutron, pion, kaon). Charged-pion production is
  comparable to Geant4's FTFP/Binary models; that is the observable for which the
  interface was validated (p+Cu / p+Pb at 3-15 GeV/c).
- Ion-ion (nucleus-nucleus) collisions are not enabled by default. UrQMD-1.3cr is the
  "cosmic ray" fork of UrQMD 1.3, trimmed for hadron+air use by shrinking its fixed
  arrays - notably the collision table (ncollmax 10000 -> 100) and the particle array
  (nmax 40000 -> 500). A dense A+A event (e.g. S+Al) overflows the 100-entry collision
  table, and the overflow bookkeeping corrupts and aborts inside UrQMD (an
  "anndec: no final state" FORTRAN stop, or a segfault). This is a capacity limit, not
  broken physics: enlarging ncollmax (in urqmd1_3/urqmd-1.3cr/colltab.f, with the
  matching size in include/G4UrQMD1_3Interface.hh) removes the crash. However,
  nucleus-nucleus was never validated for this trimmed fork.
- No coalescence / residual-nucleus de-excitation is applied, so all target nucleons
  are emitted as free particles; nucleon multiplicities are therefore higher than the
  reference Geant4 models (charged-pion production is unaffected).
- UrQMD is registered from 0 MeV, but that is a software dispatch setting, not a
  validated range. In CORSIKA it is the "low-energy" generator only relative to the
  air-shower high-energy models (below the ~80 GeV hand-off) - not a claim about the
  sub-GeV nuclear regime. Results below a few hundred MeV/nucleon need dedicated
  validation; at tens of MeV/nucleon the missing residual-nucleus / evaporation /
  coalescence physics (see above) is a real limitation.
- Not thread-safe. The UrQMD FORTRAN keeps all its state in global COMMON blocks (and is
  compiled with -fno-automatic, so even local variables are static), so it must run
  single-threaded. This example uses a serial run manager for that reason; do not use the
  UrQMD physics list in a multithreaded (MT) application.

## ACTIVATION OF CRMC INTERFACE                        

The CRMC (Cosmic Ray Monte Carlo) interface is NOT provided with Geant4 code-base.
A modified version of the CRMC interface for Geant4 applications has been kindly
prepared by Tanguy Pierog (IKP) and Andrii Tykhonov (Universite' de Geneve)
and can be obtained here:
 https://gitlab.ikp.kit.edu/AirShowerPhysics/crmc/-/tree/svn/geant4

Assuming that this special version of CRMC is installed in the subdirectory
crmc-svn-geant4/ , you need first to build it : please look at the README and
README_GEANT4_CRMC_INTERFACE files for detailed instructions on how to build it.
In short:

1. Install BOOST
2. Install HepMC (and define the corresponding environmental variable HEP_ROOT)
3. Install FASTJET (and define the corresponding environmental variable
                    FASTJET_ROOT_DIR)
4. Set the LD_LIBRARY_PATH as follows:
   ```
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HEP_ROOT}/lib:${FASTJET_ROOT_DIR}/lib
   ```
5. Source the Geant4 script geant4make.sh , e.g.
   ```
   source /your-geant4-installation-dir/share/Geant4-10.7.1/geant4make/geant4make.sh
   ```
6. Compile:
   ```
   cd crmc-svn-geant4/
   mkdir Build/ ; cd Build/   # Subdirectory where to build and install CRMC
   cmake ../
   make
   make install   # Yes, you need also to install it (in the same directory)!
   ```

After you have built CRMC you can build the Hadr02 application that uses it as follows:

1. Define the following environmental variable (in addition to the environmental
   variables defined above, needed to build CRMC):
   ```
   export G4_USE_CRMC=1
   export CRMCROOT=/your-crmc-installation-dir/crmc-svn-geant4/
   export CPATH=${CPATH}:${CRMCROOT}/Build/src:${CRMCROOT}/src
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CRMCROOT}/Build/lib
   export CRMC_CONFIG_FILE=${CRMCROOT}/Build/crmc.param
   ```
2. Compile
   ```
   cd /your-geant4/examples/extended/hadronic/Hadr02
   mkdir Build/ ; cd Build/   #  Subdirectory where to build Hadr02
   cmake -DG4_USE_CRMC=ON -DGeant4_DIR=/your-geant4-installation-dir/ ../
   make
   ```

To run the application:

1. Define the following environmental variable (besides the previous ones):
   ``` 
   export PHYSLIST=CRMC_FTFP_BERT
   ```
2. Run application:
   ```
   cd /your-geant4/examples/extended/hadronic/Hadr02/Build
   ./Hadr02 crmc.in
   ```

which runs the special "CRMC_FTFP_BERT" physics list, defined in this example,
which consists of using the standard FTFP_BERT physics list for hadrons of
kinetic energies below 100 GeV, while using CRMC above 110 GeV : in the interval
between 100 and 110 GeV, there is the transition between FTFP and CRMC (which
means that one of these two models is randomly chosen for each interaction,
with a probability which is 100% (0%) for FTFP (CRMC) at 100 GeV, and
decreases (grows) linearly to 0% (100%) for FTFP (CRMC) at 110 GeV.
Which of the MC generators of CRMC is actually used is specified in the file: \n
  `include/G4CRMCModel.hh` \n
(search for string "***LOOKHERE***" : these are the available choices:
 EPOS LHC (0) - the default - , EPOS 1.99 (1), SIBYLL 2.3c (6), and
 DPMJET 3 (12) ).

Notice that we use CRMC only for inelastic final-state of pion- , kaon- ,
proton- , neutron- and ion-nuclear interactions, whereas for the rest
(i.e. elastic and inelastic cross sections, elastic final-state interactions,
hyperon- , antihyperon- , antinucleon- and light anti-ion nuclear interactions)
we use Geant4 FTFP_BERT.

## GEOMETRY

The Target volume is a cylinder placed inside Check cylindrical volume. The 
Check volume is placed inside the World volume. The radius and the length of
the Check volume are 1 mm larger than the radius and the length of the Target.
The material of the Check volume is the same as the World material. The World
volume has the sizes 10 mm larger than that of the Target volume. Any material
from the Geant4 database can be defined. The default World  material is
G4Galactic and the default  Target material is aluminum. The Target is
subdivided into a number of equal slices. The following UI commands are available to
modify the geometry:

```
/testhadr/TargetMat     G4_Pb
/testhadr/WorldMat      G4_AIR
/testhadr/TargetRadius  10 mm
/testhadr/TargetLength  20 cm
/testhadr/NumberDivZ    200
```

Beam direction coincides with the target axis and is Z axis in the global
coordinate system. G4ParticleGun is used as a primary generator. The energy 
and the type of the beam can be defined via standard UI commands

```
/gun/energy   150 GeV
/gun/particle ion
/gun/ion 6 12
```

Default beam position is -(targetHalfLength + 5*mm) and direction along Z axis.
Beam position and direction can be changed by gun UI commands:

```
/gun/position  1 10 3 mm
/gun/direction 1 0 0
```

however, position command is active only if before it the flag is set

```
/testhadr/DefaultBeamPosition false
```
 
## SCORING

The scoring is performed with the help of UserStackingAction class and a
sensitive detector class associated with a target slice. 
Each secondary particle is scored by the StackingAction.  In
the StackingAction it is also possible to kill all or only EM (e+, e-, gamma)
secondary particles 

```
/testhadr/killAll  
/testhadr/KillEM
```

To control running the following options are available:

```
/run/printProgress 10
```

## PHYSICS

PhysicsList of the application uses components, which are distributed with
Geant4 in /geant4/physics_lists subdirectory. 

Reference Physics Lists are used and the environment variable PHYSLIST should 
be defined. 

Additionally it is possible to add ion-ion interactions using UI command

```
/testhadr/ionPhysics   HIJING
/testhadr/ionPhysics   UrQMD
```

## VISUALIZATION

The vis.mac file can be used as an example of visualization.

## HISTOGRAMS

It is possible to choose the format of the output file with 
histograms using UI command:

```
/testhadr/histo/fileName   name
/testhadr/histo/fileType   type
```

The following types are available: root, xml(aida). They will be
stored in the file "name.root", or "name.xml".


All histograms are normalized to the number of events.

