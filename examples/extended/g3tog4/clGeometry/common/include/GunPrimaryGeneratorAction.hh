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
/// \file GunPrimaryGeneratorAction.hh
/// \brief Definition of the GunPrimaryGeneratorAction class

#ifndef GunPrimaryGeneratorAction_h
#define GunPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class G4ParticleGun;
class G4Event;

/// \ingroup primary_generator
/// \brief The primary generator class with particle gun
///
/// \author I. Hrivnacova; IPN Orsay

class GunPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    GunPrimaryGeneratorAction(
      const G4String& particleName = "geantino",
      G4double energy = 1.*CLHEP::MeV,
      G4ThreeVector position= G4ThreeVector(0,0,0),
      G4ThreeVector momentumDirection = G4ThreeVector(0,0,1));    
    ~GunPrimaryGeneratorAction();

    // methods
    virtual void GeneratePrimaries(G4Event*);

  private:
    // data members
    G4ParticleGun*  fParticleGun;
};

#endif
