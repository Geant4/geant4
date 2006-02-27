// $Id: QPrimaryGeneratorAction.hh,v 1.1 2006-02-27 10:05:24 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QPrimaryGeneratorAction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef Q_PRIMARY_GENERATOR_ACTION_H
#define Q_PRIMARY_GENERATOR_ACTION_H

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;

// ====================================================================
//
// class definition
//
// ====================================================================
class QPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
  // use G4 particle gun
  G4ParticleGun* particleGun;

public:
  QPrimaryGeneratorAction();
  ~QPrimaryGeneratorAction();

  G4ParticleGun* GetParticleGun() const;

  virtual void GeneratePrimaries(G4Event* anEvent);
};

// ====================================================================
//   inline functions
// ====================================================================
inline G4ParticleGun* QPrimaryGeneratorAction::GetParticleGun() const
{  return particleGun; }

#endif
