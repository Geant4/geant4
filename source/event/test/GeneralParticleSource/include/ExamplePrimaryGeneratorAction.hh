#ifndef ExamplePrimaryGeneratorAction_h
#define ExamplePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4GeneralParticleSource;
class G4Event;

class ExamplePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExamplePrimaryGeneratorAction();
    ~ExamplePrimaryGeneratorAction();

  // DEBUG arrays
  G4int debugx[100], debugy[100], debugz[100];
  G4int debugpx[100], debugpy[100], debugpz[100];
  G4int debugtheta[100], debugphi[100];
  G4int debugenergy[100];
  G4double debug_energy_step;
  G4double DebugXmin, DebugXmax, DebugYmin, DebugYmax;
  G4double DebugZmin, DebugZmax, DebugXStep, DebugYStep;
  G4double DebugZStep;

  G4String SourceType;
  G4String SourceShape;
  G4double radius;
  G4double halfx;
  G4double halfy;
  G4double halfz;
  G4ThreeVector centre;

    // for use with energy
  G4String EneDisType ;
  G4String InterpolationType;
  G4double emin;
  G4double emax;

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4GeneralParticleSource* particleGun;
  
};

#endif



