//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: QPrimaryGeneratorAction.hh,v 1.2 2006-06-04 21:35:59 kmura Exp $
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
