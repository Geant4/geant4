\page ExampleAuNP Example AuNP

\author D. Sakata (a)                             \n
(a) email: sakata.dousatsu@qst.go.jp, dosatsu.sakata@cern.ch \n

## INTRODUCTION

The AuNP example simulates the track-structure of electrons in microscopic gold volume.
The example also simulates that in liquid water medium surrounding the gold volume.

This example is provided by the Geant4-DNA collaboration.

These processes and models are further described at:\n
http://geant4-dna.org                               \n\n

Any report or published results obtained using the Geant4-DNA software shall
cite the following Geant4-DNA collaboration publications:\n
Phys. Med. 31 (2015) 861-874                             \n
Med. Phys. 37 (2010) 4692-4708                           \n\n

We also suggest these other references related to this example:\n
Phys. Med. 63 (2019), 98-104                                   \n
Med. Phys. 45 (2018), 2230-2242                                \n
J. App. Phys. 120 (2016), 244901                               \n\n

The AuNP example simulates the track-structure of electrons
in microscopic gold volume.The example also simulates that
in liquid water medium surrunding the gold volume.

The details of the physics models are described in the following paper:
J. Appl. Phys. 120 (2016) 244901

The details of the application are described in the following paper:
Phys. Med. 63 (2019) 98-104
Med. Phys. 45(5) (2018) 2230-2242

## GEOMETRY SET-UP

The geometry is a sphere (World) made of liquid water containing a smaller
sphererical target volume of Gold.

The default geometry is constructed in DetectorConstruction class.

## PHYSICS LIST

The default physics list is constructed in PhysicsList class.

## DETECTOR RESPONSE : Scorers

This scorer computes following quantities.

- the energy spectra of secondary particles generated in AuNP
- the energy spectra of secondary particles at AuNP surface
- the energy spectra of secondary particles generated in liquid water
- the energy deposit and the position in the absorber surrunding AuNP

Run::RecordEvent(), called at end of event, collects informations
event per event from the hits collections, and accumulates statistic for
RunAction::EndOfRunAction().

In multi-threading mode the statistics accumulated per workers is merged
to the master in Run::Merge().

The information is scored in a ROOT ntuple file AuNP.root.
  
## HOW TO RUN THE EXAMPLE

This example shows:
- how to use the Geant4-DNA processes for gold
- how to count and save occurrences of processes

The code can be compiled with cmake.
It works in MT mode.

Two user macro files can be used:
```
./AuNP AuNP.mac
```

## SIMULATION OUTPUT AND RESULT ANALYSIS

The output results consists in an AuNP.root file, containing for the run:
- the energy spectra of secondary particles generated in AuNP
- the energy spectra of secondary particles at AuNP surface
- the energy spectra of secondary particles generated in liquid water
- the energy deposit and the position in the absorber surrounding AuNP

This file can be easily analyzed using for example the provided ROOT macro
file plot.C. The plot.C is for for study of physics only. 
Make sure the chemistry is disactivated as "/chem/activate False".

to do so :
- be sure to have ROOT installed on your machine
- be sure to be in the build directory of the AuNP example
- launch ROOT by typing root
- under your ROOT session, type in: .X plot.C to execute the macro file
- alternatively you can type directly under your session: root plot.C
