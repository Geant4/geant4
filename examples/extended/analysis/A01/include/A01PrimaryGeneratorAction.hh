// $Id: A01PrimaryGeneratorAction.hh,v 1.1 2002-11-13 07:18:53 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef A01PrimaryGeneratorAction_h
#define A01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4ParticleDefinition;
class A01PrimaryGeneratorMessenger;

class A01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    A01PrimaryGeneratorAction();
    virtual ~A01PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* particleGun;
    A01PrimaryGeneratorMessenger* gunMessenger;
    G4ParticleDefinition* positron;
    G4ParticleDefinition* muon;
    G4ParticleDefinition* pion;
    G4ParticleDefinition* kaon;
    G4ParticleDefinition* proton;
    G4double momentum;
    G4double sigmaMomentum;
    G4double sigmaAngle;
    G4bool randomizePrimary;

  public:
    inline void SetMomentum(G4double val) { momentum = val; }
    inline G4double GetMomentum() const { return momentum; }
    inline void SetSigmaMomentum(G4double val) { sigmaMomentum = val; }
    inline G4double GetSigmaMomentum() const { return sigmaMomentum; }
    inline void SetSigmaAngle(G4double val) { sigmaAngle = val; }
    inline G4double GetSigmaAngle() const { return sigmaAngle; }
    inline void SetRandomize(G4bool val) { randomizePrimary = val; }
    inline G4bool GetRandomize() const { return randomizePrimary; }
};

#endif


