// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: gogdmlPrimaryGeneratorAction.hh,v 1.1.1.1 2002-05-31 00:34:43 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef gogdmlPrimaryGeneratorAction_h
#define gogdmlPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class gogdmlPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    gogdmlPrimaryGeneratorAction();
    ~gogdmlPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


