// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst18PrimaryGeneratorAction.hh,v 1.4 2000-06-14 12:33:12 flei Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst18PrimaryGeneratorAction_h
#define Tst18PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4IonTable.hh"
#include "RadioactiveDecayGun.hh"

//class G4ParticleGun;

class G4Event;

class Tst18PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst18PrimaryGeneratorAction();
    ~Tst18PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  protected:

     RadioactiveDecayGun *theParticleGun;

};

#endif


