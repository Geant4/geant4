#ifndef hTestPrimaryGeneratorAction_h
#define hTestPrimaryGeneratorAction_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestPrimaryGeneratorAction
//  
// Description: Generate primary beam 
//
// Authors:    0.6.04.01 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestEventAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestPrimaryGeneratorMessenger;

class hTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

  // The constructor defines a ParticleGun object, which allows 
  // shooting a beam of particles through the experimental set-up.
    hTestPrimaryGeneratorAction(hTestEventAction*);

  //The destructor. It deletes the ParticleGun.
    ~hTestPrimaryGeneratorAction();

  public:

  //It reads the parameters of the primary particles.
  //Generates the primary event via the ParticleGun method.
    void GeneratePrimaries(G4Event* anEvent);

  //Get/Set methods
    void SetBeamX(G4double val) {x0 = val;};
    void SetBeamY(G4double val) {y0 = val;};
    void SetBeamZ(G4double val) {z0 = val;};
    void SetBeamSigmaX(G4double val) {sigmaX = val;};
    void SetBeamSigmaY(G4double val) {sigmaY = val;};
    void SetBeamSigmaZ(G4double val) {sigmaY = val;};
    void SetBeamSigmaE(G4double val) {sigmaE = val;};
    void SetBeamMinCosTheta(G4double val) {minCosTheta = val;};
    void SetVerbose(G4int val) {verbose = val;};
    G4ThreeVector GetBeamPosition() const {return position;}; 
    G4ThreeVector GetBeamDirection() const {return direction;};
    G4ThreeVector GetBeamEnergy() const {return energy;};

  private:

    void InitializeMe();

    hTestEventAction* theEvent;
    G4ParticleGun* particleGun;
    hTestPrimaryGeneratorMessenger* theMessenger;

    G4int counter;
    G4int verbose;
    G4double x0, y0, z0;
    G4double sigmaX, sigmaY, sigmaZ;
    G4double sigmaE;
    G4double energy;
    G4double minCosTheta;
    G4ThreeVector position;
    G4ThreeVector direction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif






