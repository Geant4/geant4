// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14PrimaryGeneratorAction.hh,v 1.2 1999-06-14 14:28:34 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst14PrimaryGeneratorAction_h
#define Tst14PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst14PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst14PrimaryGeneratorAction();
    ~Tst14PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


