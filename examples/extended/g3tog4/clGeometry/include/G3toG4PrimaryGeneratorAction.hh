// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4PrimaryGeneratorAction.hh,v 1.1 2000-07-24 11:23:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G3toG4PrimaryGeneratorAction_h
#define G3toG4PrimaryGeneratorAction_h 1

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"

class G3toG4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    G3toG4PrimaryGeneratorAction();
    ~G3toG4PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
    G4ThreeVector GetRandomDirection();
};

#endif


