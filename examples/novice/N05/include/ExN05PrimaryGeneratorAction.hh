// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05PrimaryGeneratorAction.hh,v 1.2 1999-12-15 14:49:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN05PrimaryGeneratorAction_h
#define ExN05PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class ExN05PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExN05PrimaryGeneratorAction();
    ~ExN05PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    G4ParticleGun* GetParticleGun();

  private:
    G4ParticleGun* particleGun;
};

#endif


