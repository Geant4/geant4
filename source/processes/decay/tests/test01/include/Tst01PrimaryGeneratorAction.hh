// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst01PrimaryGeneratorAction.hh,v 1.1 2001-02-08 08:41:45 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst01PrimaryGeneratorAction_h
#define Tst01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class Tst01ParticleGun;
class G4Event;

class Tst01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst01PrimaryGeneratorAction();
    ~Tst01PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    Tst01ParticleGun* particleGun;
};

#endif


