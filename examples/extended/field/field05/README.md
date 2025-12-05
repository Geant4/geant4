\page Examplefield05 Example field05

 This example checks so-called "spin-frozen" condition
 There is a good example article hep-ph/0012087v1.
 This article discusses about how to cancel the muon g-2 precession by
 applying an electric field.

 - 1) beta is muon velocity,
 - 2) B is an uniform magnetic field and vec{beta}.vec{B}=0,
    "." means scalar product,
 - 3) Radial electric field (E) in the lab frame and vec{beta}.vec{E}=0,
 - 4) a=(g-2)/2 is muon anomalous magnetic moment.

 The required electric field to cancel the g-2 precession is,
```cpp
E = a*B*light_c*gamma**2*beta
```

 In case of gamma=5 and B=0.24 Tesla, the required electric field is
```cpp
E = 2 MV/m
```

 "Spin-frozen" happens when spin rotation cycle and muon rotation cycle
 are same. In this case, both cycles should be 149.5 nsec.

 See also:
 http://research.kek.jp/people/hiromi/MyHomePage/G-2_work_files/SpinStudyinEMfieldByGeant4.pdf

 Credit goes to Hiromi Iinuma from KEK.

## main()

 See main() in field05.cc.

## GEOMETRY DEFINITION

 As simple world G4Box with a G4ElectroMagneticField      \n
 propagating both spin and momentum (G4EqEMFieldWithSpin) \n
 with G4ClassicalRK4(fEquation,12) and                    \n
 Bz = 0.24*tesla;                                         \n
 Er = 2.113987E+6*volt/m;

## AN EVENT: THE PRIMARY GENERATOR

 Use mu+ G4ParticleGun with Pmu = 517.6*MeV/c            \n
 and aligned spin and momentum direction

## PHYSICS

```cpp
RegisterPhysics(new G4SpinDecayPhysics());
RegisterPhysics(new G4StepLimiterPhysics());
```

 - G4SpinDecayPhysics defines muon decay modes with spin,
 - G4StepLimiterPhysics defines G4StepLimiter and G4UserSpecialCuts.

## User Action Classes

- F05SteppingAction:

  G4Exception when the cosine of the angle between
  the spin and the momentum is < (1.-1.E-7)

## HOW TO START ?

 - Execute field05 in 'batch' mode from macro files e.g.
```
% ./field05 field05.in > field.out &
```

 - Execute field05 in 'interactive' mode with visualization e.g.
```
% ./field05
....
Idle> type your commands, for example:
Idle> run/beamOn 1
....
```