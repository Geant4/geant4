// $Id: ParticleGunAction.hh,v 1.2 2006-04-25 10:29:52 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   ParticleGunAction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef PARTICLE_GUN_ACTION_H
#define PARTICLE_GUN_ACTION_H

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;

// ====================================================================
//
// class definition
//
// ====================================================================
class ParticleGunAction : public G4VUserPrimaryGeneratorAction {
private:
  // use G4 particle gun
  G4ParticleGun* particleGun;

public:
  ParticleGunAction();
  ~ParticleGunAction();

  G4ParticleGun* GetParticleGun() const;

  virtual void GeneratePrimaries(G4Event* anEvent);
};

// ====================================================================
//   inline functions
// ====================================================================
inline G4ParticleGun* ParticleGunAction::GetParticleGun() const
{  return particleGun; }

#endif
