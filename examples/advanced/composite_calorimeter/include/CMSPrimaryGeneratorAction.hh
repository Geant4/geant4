///////////////////////////////////////////////////////////////////////////////
// File: CMSPrimaryGeneratorAction.hh
// Date: 06/98
// Description: Sets up particle beam
//
//     By default 1 pi+  is shot from (0,0,0)
//     in (1,1,0.1) direction at 100 GeV
//     Use /gun/... commands to modify energy,origin,direction at run time.
//     or/and
//         /OSCAR/generator/random true/false to have random direction
//         /OSCAR/generator/scan   true/false to scan in eta/phi
//     Use 
//         /OSCAR/generator/minEnergy
//         /OSCAR/generator/maxEnergy
//         /OSCAR/generator/minPhi 
//         /OSCAR/generator/maxPhi 
//         /OSCAR/generator/minEta
//         /OSCAR/generator/maxEta 
//     to set the range in energy and direction of particles shot at random.
//     Use 
//         /OSCAR/generator/stepsPhi
//         /OSCAR/generator/stepsEta
//     to set number of steps in Phi and Eta for the scan
//
//
// Last modified: 08/98 I.G. (updated to beta version)
//                06/08/99 V.Lefebure -> add possibility to use pythia events
//                08/09/99 I.G.       -> Add gunMessenger. Pythia file clean up
//                12/10/99 V.L.       -> add comments
//		  14/10/99 V.L.	      -> add solid angle restriction
//		  15/10/99 V.L.       -> add gun messenger options
//                18/04/00 P.A., S.B. -> Extended functionality
//                   10/01 P.Arce use COBRA GeneratorInterface
///////////////////////////////////////////////////////////////////////////////


#ifndef CMSPrimaryGeneratorAction_h
#define CMSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4VPrimaryGenerator.hh"

#include "CMSPrimaryGeneratorMessenger.hh"

enum generatorInputType {singleFixed, singleRandom, singleScan};

class CMSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  CMSPrimaryGeneratorAction();
  ~CMSPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  
public:
  void SetVerboseLevel(G4int val);
  void SetRandom(G4String val);
  void SetScan(G4String val);
  void SetMinimumEnergy(G4double p);
  void SetMaximumEnergy(G4double p);
  void SetMinimumPhi(G4double p);
  void SetMaximumPhi(G4double p);
  void SetStepsPhi(G4int val);
  void SetMinimumEta(G4double p);
  void SetMaximumEta(G4double p);
  void SetStepsEta(G4int val);
  void SetGunPosition(const G4ThreeVector & pos) const;
  void SetRunNo(G4int val);

public:    
  G4ThreeVector GetParticlePosition() {return particleGun->GetParticlePosition();}
  G4double GetParticleEnergy() {return particleGun->GetParticleEnergy();}

private:
  CMSPrimaryGeneratorMessenger* gunMessenger;
  G4ParticleGun* particleGun;
  generatorInputType generatorInput;

  G4int verboseLevel;
  G4int n_particle;
  G4String particleName;
  G4double particleEnergy;
  G4ThreeVector particlePosition;
  G4ThreeVector particleDir;

  G4double energyMin,energyMax;
  G4double etaMin,etaMax;
  G4double phiMin,phiMax;
  G4int etaSteps,phiSteps;

  G4int isInitialized;
  G4double etaValue, phiValue;
  G4int scanSteps;

private:
  void initialize();
  void print(G4int val);

};

#endif



