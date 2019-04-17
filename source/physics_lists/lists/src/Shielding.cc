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
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2010 Tatsumi Koi, Gunter Folger
//
//   created from FTFP_BERT
//
// Modified:
// 05.08.2014 K.L.Genser: added provision for Hadronic Physics Variant M
// 16.08.2010 H.Kurashige: Remove inclusion of G4ParticleWithCuts 
// 26.04.2011 T.Koi: Add G4RadioactiveDecayPhysics
// 16.10.2012 A.Ribon: Use new default stopping
// 07.11.2013 T.Koi: Add IonElasticPhysics, Set proton cut to 0 to generate
//                   low energy recoils and activate production of fission
//                   fragments 
//
//----------------------------------------------------------------------------
//

#include "Shielding.hh"
#include "globals.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmParameters.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonQMDPhysics.hh"
#include "G4IonElasticPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsLEND.hh"
#include "G4ParticleHPManager.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsShielding.hh"
#include "G4HadronPhysicsShieldingLEND.hh"
#include <CLHEP/Units/SystemOfUnits.h>

Shielding::Shielding(G4int verbose, const G4String& n_model, 
                     const G4String& HadrPhysVariant )
{
  G4String LEN_model = n_model; 
  size_t find = LEN_model.find("LEND__");
  G4String evaluation;
  if ( find != G4String::npos )
  {
     evaluation=LEN_model;
     evaluation.erase(0,find+6);
     LEN_model="LEND";
  }

  G4DataQuestionaire it(photon, neutron, radioactive);
  G4cout << "<<< Geant4 Physics List simulation engine: Shielding"
         << HadrPhysVariant << G4endl;
  if ( LEN_model=="LEND" ) G4cout << 
    "<<< LEND will be used for low energy neutron and gamma projectiles"
    << G4endl;

  defaultCutValue = 0.7*CLHEP::mm;  
  SetCutValue(0, "proton");  
  SetVerboseLevel(verbose);

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(verbose));

  // Synchroton Radiation & GN Physics
  G4EmExtraPhysics* emExtraPhysics = new G4EmExtraPhysics(verbose);
  if ( LEN_model == "LEND" ) {
     // Use LEND model for Gamma Nuclear 
     emExtraPhysics->LENDGammaNuclear(true);
  }
  RegisterPhysics( emExtraPhysics );

  // Decays 
  RegisterPhysics( new G4DecayPhysics(verbose) );
  RegisterPhysics( new G4RadioactiveDecayPhysics(verbose) );

  // Hadron Elastic scattering
  if ( LEN_model == "HP" ) 
  {
     RegisterPhysics( new G4HadronElasticPhysicsHP(verbose) );
  }
  else if ( LEN_model == "LEND" ) 
  {
     RegisterPhysics( new G4HadronElasticPhysicsLEND(verbose,evaluation));
     G4DataQuestionaire itt(lend);
  }
  else 
  {
     G4cout << "Shielding Physics List: Warning!" <<G4endl;
     G4cout << "\"" << LEN_model 
            << "\" is not valid for the low energy neutorn model." <<G4endl;
     G4cout << "Neutron HP package will be used." <<G4endl;
     RegisterPhysics( new G4HadronElasticPhysicsHP(verbose) );
  } 

  G4VPhysicsConstructor* hpc;
  // Hadron Physics HP or LEND 
  if (HadrPhysVariant == "M") {
     hpc = new G4HadronPhysicsShielding("hInelastic Shielding", verbose, 
                                        9.5*CLHEP::GeV, 9.9*CLHEP::GeV);
  } else {
     hpc = new G4HadronPhysicsShielding("hInelastic Shielding", verbose, 
                                        4.0*CLHEP::GeV, 5.0*CLHEP::GeV);
  }

  if ( LEN_model == "LEND" ) {
     delete hpc;
     if (HadrPhysVariant == "M") {
        hpc = new G4HadronPhysicsShieldingLEND("hInelastic ShieldingLEND", 
                                 verbose, 9.5*CLHEP::GeV, 9.9*CLHEP::GeV);
     } else {
        hpc = new G4HadronPhysicsShieldingLEND("hInelastic ShieldingLEND", 
                                 verbose, 4.0*CLHEP::GeV, 5.0*CLHEP::GeV);
     }
  } else {
     //G4cout << "Shielding Physics List: Warning." <<G4endl;
     //G4cout << "Name of Low Energy Neutron model " << LEN_model 
     //         << " is invalid." <<G4endl;
     //G4cout << "Will use neutron HP package." <<G4endl;
  }
  RegisterPhysics( hpc );

  if ( LEN_model == "HP" ) {
     //Activate prodcuton of fission fragments in neutronHP
     G4ParticleHPManager::GetInstance()->SetProduceFissionFragments( true );
  }

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(verbose) );

  // Ion Physics
  RegisterPhysics( new G4IonElasticPhysics(verbose));
  RegisterPhysics( new G4IonQMDPhysics(verbose));  
}

Shielding::~Shielding()
{}
