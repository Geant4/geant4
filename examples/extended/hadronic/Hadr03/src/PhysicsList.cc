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
/// \file hadronic/Hadr03/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
<<<<<<< HEAD
// $Id: PhysicsList.cc 86125 2014-11-07 11:07:29Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4HadronElasticPhysicsHP.hh"
<<<<<<< HEAD
=======
#include "G4HadronElasticPhysicsXS.hh"

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronPhysicsShielding.hh"

#include "G4IonElasticPhysics.hh"
<<<<<<< HEAD
#include "G4IonPhysics.hh"
=======
#include "G4IonPhysicsXS.hh"
#include "G4IonQMDPhysics.hh"
#include "G4IonPhysicsPHP.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
#include "G4IonINCLXXPhysics.hh"

#include "GammaNuclearPhysics.hh"
#include "GammaNuclearPhysicsLEND.hh"

#include "G4RadioactiveDecayPhysics.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
:G4VModularPhysicsList()
{
  G4int verb = 1;  
  SetVerboseLevel(verb);
  
  //add new units for cross sections
  // 
  new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);  
  
  // Hadron Elastic scattering
  //
  RegisterPhysics( new G4HadronElasticPhysicsHP(verb));
  ///RegisterPhysics( new G4HadronElasticPhysicsXS(verb));  

  // Hadron Inelastic physics
  //
<<<<<<< HEAD
  RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));
  ////RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));
  ////RegisterPhysics( new G4HadronInelasticQBBC(verb));        
=======
  ////RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));
  RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));
  ////RegisterPhysics( new G4HadronPhysicsQGSP_BIC_AllHP(verb));
  ////RegisterPhysics( new G4HadronPhysicsQGSP_BIC(verb));  
  ////RegisterPhysics( new G4HadronInelasticQBBC(verb));
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));
  ////RegisterPhysics( new G4HadronPhysicsShielding(verb));
    
  // Ion Elastic scattering
  //
  RegisterPhysics( new G4IonElasticPhysics(verb));  
  
  // Ion Inelastic physics
  //
<<<<<<< HEAD
  RegisterPhysics( new G4IonPhysics(verb));
=======
  RegisterPhysics( new G4IonPhysicsXS(verb));
  ////RegisterPhysics( new G4IonPhysicsPHP(verb));
  ////RegisterPhysics( new G4IonQMDPhysics(verb));
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  ////RegisterPhysics( new G4IonINCLXXPhysics(verb));
    
  // Gamma physics
  //
<<<<<<< HEAD
  RegisterPhysics( new GammaPhysics("gamma"));
=======
  RegisterPhysics( new GammaNuclearPhysics("gamma"));
  ////RegisterPhysics( new GammaNuclearPhysicsLEND("gamma")); 
  
  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
<<<<<<< HEAD
  SetCutValue(0*mm, "proton");
=======
   SetCutValue(0.*mm, "proton");
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
