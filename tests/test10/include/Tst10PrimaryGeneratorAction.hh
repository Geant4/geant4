// $Id: Tst10PrimaryGeneratorAction.hh,v 1.1 1999-01-08 16:35:32 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPrimaryGeneratorAction
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Tst10PrimaryGeneratorAction_h
#define Tst10PrimaryGeneratorAction_h 1

#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst10PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst10PrimaryGeneratorAction();
    ~Tst10PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
  private:
    G4ParticleGun* particleGun;
		G4ThreeVector GetRandomPolarization( G4ThreeVector direction );
		G4ThreeVector GetRandomDirection( );
		G4ThreeVector GetRandomPosition( );
};

#endif


