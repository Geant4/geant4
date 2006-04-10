// -------------------------------------------------------------------
// $Id: MicrobeamPrimaryGeneratorAction.hh,v 1.2 2006-04-10 14:47:31 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamPrimaryGeneratorAction_h
#define MicrobeamPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

#include "MicrobeamDetectorConstruction.hh"

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MicrobeamPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  MicrobeamPrimaryGeneratorAction(MicrobeamDetectorConstruction*);    
  ~MicrobeamPrimaryGeneratorAction();
  
  void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun*           particleGun;
  MicrobeamDetectorConstruction*    Detector;     
};

#endif


