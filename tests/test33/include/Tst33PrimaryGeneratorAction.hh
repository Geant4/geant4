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
//
// $Id: Tst33PrimaryGeneratorAction.hh,v 1.2 2002-11-20 13:09:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33PrimaryGeneratorAction
//
// Class description:
//
// Create 10 MeV neutrons.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33PrimaryGeneratorAction_hh
#define Tst33PrimaryGeneratorAction_hh Tst33PrimaryGeneratorAction_hh 

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst33PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  Tst33PrimaryGeneratorAction();
  virtual ~Tst33PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* anEvent);

private:
  Tst33PrimaryGeneratorAction(const Tst33PrimaryGeneratorAction &);
  Tst33PrimaryGeneratorAction &operator=(const Tst33PrimaryGeneratorAction &);
  G4ParticleGun* fParticleGun;
};

#endif
