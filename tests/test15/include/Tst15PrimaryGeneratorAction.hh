// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst15PrimaryGeneratorAction.hh,v 1.1 1999-11-18 14:48:13 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst15PrimaryGeneratorAction_h
#define Tst15PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst15PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst15PrimaryGeneratorAction();
    ~Tst15PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


