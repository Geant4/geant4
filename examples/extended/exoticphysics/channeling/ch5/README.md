\page Examplech5 Example ch5

\author  G. Paternò, A. Sytov - INFN Ferrara Division (Italy) \n
paterno@fe.infn.it, sytov@fe.infn.it

### INTRODUCTION
Example ch5 is an application for simulating a positron source. 
Although the conventional approach based on an **amorphous target** is possible, 
the application is primarily designed to simulate positron sources based on oriented crystals. 
In the latter case, both the **single-crystal** and the **hybrid scheme** 
can be investigated [[1]](#1). 

### DESCRIPTION
One or two main volumes can be present in the setup, depending if the user wants to consider 
a conventional/single-oriented-crystal scheme [[2]](#2) or the hybrid scheme [[3]](#3). 
In the hybrid scheme, the first volume is an oriented crystal (typically along a crystalline axis)
that serves as a radiator, whereas the second volume is a randomly oriented crystal 
(equivalent to an amorphous volume) where the photons emitted by the radiator are converted 
into positrons.
The converter can be composed of small spheres (the so called _granular target_ - GT) 
so as to reduce the energy deposition and the thermomechanical stress. 
In addition, from _one to three scoring screens_ are present to score the particles
leaving or enetering the aformentioned volumes. In particular the scoring screen identified
with number 2 is positioned just downstram of the radiator, while scoring screens 0 and 1 
are positioned just upstream and downstream of the converter, respectively. 
In a conventional or a single crystal scheme, only the scoring screen 0 is present 
and it is automatically positioned just downstream of the single volume positioned.

An **advanced hybrid scheme that includes an ideal bending magnet or a collimator**
to remove the charged particle or limit the number of particles impinging on the converter, 
respectively, **can also be considered** [[1]](#1). 

Through a set of custom macro commands, the user can define the geometry ad the scoring strategy.
A description of all the available options is provided in _run.mac_ (inside the macros folder). 
As an example, the Orientetional Coherent (OC) effects (including radiation) in crystals 
(enabled by G4ChannelingFastSimModel), which by deafult are activated, can be deactivated 
through the command: `/crystal/setOCeffects false`.

The back surface of the radiator crystal is placed at z=0 (with z as the beam direction), 
while the front position of the possible converter can be set up via macro.

Various macros are available to simulate different configurations: run_conventional.mac, 
run_single_crystal.mac, run_hybrid.mac, and run_hybrid_granular_target.mac 
for a conventional, single oriented crystal and hybrid with solid or granular target, 
respectively. The parameters set in these macros come from the study carried out for 
the positron source of FCC-ee [[2]](#2). 
However, they can be changed to investigate different cases.

The output is recorded into a root file whose name can be set by macro 
(default is output/output.root), as a set of ntuples.

The ntuple "scoring_ntuple" is used to score the features of the particles impinging 
on the scoring screens. It contains the following variables (columns):

"screenID", "particle", "x", "y", "px", "py", "pz", "t", "eventID"

which represents:
- the screen ID (column 0),
- the particle name (column 1),
- the impinging x,y coordinates in mm (columns 2,3),
- the momentum components (MeV) of the particle (columns 4-6),
- the time of arrival of the particle in ns (column 7),
- the event ID (column 9).

The ntuple "edep_rad" and "edep_conv" are used to store the energy deposited (MeV) 
in radiator and converter, respectively, thus they contain simply the variables:
"edep", "eventID"

The ntuple "edep_spheres" is instead used to store the energy deposited (MeV) inside the spheres 
of a granular target/converter. It contains the variables:
"volumeID", "edep", "eventID"
where volumeID identify the single sphere inside the target. This ntuple is filled only if the 
target is indeed granular (it can be activated through the command /det/setGranularConverter true).

Finally, the ntuple "scoring_ntuple2" is used to score the features of the particles leaving 
the radiator or the target/converter. It contains the following variables (columns):

"particle", "x", "y", "z", "px", "py", "pz", "t", "eventID", "trackID"

which represents:
- the particle name (column 0),
- the impinging x,y,z coordinates in mm (columns 1-3),
- the momentum components (MeV) of the particle (columns 4-6),
- the time of arrival of the particle in ns (column 7),
- the event ID (column 8),
- the track ID (column 9).

The three-dimensional distributions of energy deposition in the converter 
(radiator if the converter is not present) can be scored through the standard 
box mesh scorer defined in the attached macros.

To visualize these data one should use the python notebook analysis_ch5.ipynb.

Once the example is build, an interactive session with the graphic user intergace (GUI)
showing a defualt geometry set through macro geom.mac can be run by simply typing
`./ch5` in a terminal after moving inside the build directory.

### REFERENCES
<a id="1">[1]</a> M. Soldani, et al. NIM A 1058 (2024): 168828 (https://doi.org/10.1016/j.nima.2023.168828).

<a id="2">[2]</a>  F. Alharthi et al. NIM A 1075 (2025): 170412 (https://doi.org/10.1016/j.nima.2025.170412).

<a id="3">[3]</a>  N. Canale et al. NIM A 1075 (2025): 170342 (https://doi.org/10.1016/j.nima.2025.170342). 

