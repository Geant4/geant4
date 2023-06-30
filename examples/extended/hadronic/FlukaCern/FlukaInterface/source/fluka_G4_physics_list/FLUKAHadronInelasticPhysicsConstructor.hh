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
// Construct hadron inelastic physics processes with FLUKA.CERN XS and models.
//
// NB 1: Neutron HP treatment is the one from G4.
//
// NB 2: Neutron capture and fission processes are also defined here in this constructor 
// (as is done in the G4 physics lists), and are from G4.
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef FLUKA_HADRON_INELASTIC_PHYSICS_CONSTRUCTOR_HH
#define FLUKA_HADRON_INELASTIC_PHYSICS_CONSTRUCTOR_HH


// G4
#include "CLHEP/Units/SystemOfUnits.h"
#include "globals.hh"
#include "G4VPhysicsConstructor.hh"


class FLUKAHadronInelasticPhysicsConstructor final : public G4VPhysicsConstructor {

public:
  FLUKAHadronInelasticPhysicsConstructor(G4int verbose = 1);

  virtual void ConstructParticle() override;
  virtual void ConstructProcess() override;


private:
  G4double fNeutronHPMaxE = 20*CLHEP::MeV;
  G4double fNeutronFLUKAMinE = 20*CLHEP::MeV;
};


#endif
#endif // G4_USE_FLUKA
