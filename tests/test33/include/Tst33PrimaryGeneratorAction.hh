//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Tst33PrimaryGeneratorAction.hh,v 1.3 2006-06-29 21:59:47 gunter Exp $
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
