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
// $Id: A01PhysicsList.cc,v 1.8 2006-06-29 16:33:03 gunter Exp $
// --------------------------------------------------------------
//
// 28-Jan-04 Add QGSP_BERT and QGSP_BIC for hadronic lists. T. Koi
// 22-Nov-04 Comment out QGSP_BERT and QGSP_BIC
//           Output Notificaiton message             
//           All Particles are created in GeneralPhysics 

#include "A01PhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>

#include "A01GeneralPhysics.hh"
#include "A01EMPhysics.hh"
#include "A01MuonPhysics.hh"
#include "A01HadronPhysics.hh"
#include "A01IonPhysics.hh"

//#include "HadronPhysicsQGSP_BERT.hh"
//#include "HadronPhysicsQGSP_BIC.hh"

A01PhysicsList::A01PhysicsList():  G4VModularPhysicsList()
{

  G4cout << "You are using the A01PhysicsList" << G4endl;
  G4cout << "Full set of particles (barions bosons and mesons) will be created and" << G4endl;
  G4cout << "Standard EM Physics and Low & High Energy parameterized models will be applied." << G4endl;
  G4cout << "A01PhysicsList is optimized for robustness" << G4endl;
  G4cout << "and not for any particular usage." << G4endl;
  G4cout << "For the hadronic physics, educated guesses of physics list are prepared for various use cases." << G4endl;
  G4cout << "When you will start REAL calculations for your own interest," << G4endl;
  G4cout << "please consider the usage of hadronic_lists instead of A01PhysicsLists." << G4endl;
  G4cout << "More information can also be found from the Geant4 HyperNews." << G4endl;
  G4cout << "http://geant4-hn.slac.stanford.edu:5090/Geant4-HyperNews/index" << G4endl;
  G4cout << "" << G4endl;

  // default cut value  (1.0mm)
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics ( Create ALL Particle and apply Decay )
  RegisterPhysics( new A01GeneralPhysics("general") );

  // EM Physics ( Apply related Processes to gamma and e-/+)
  RegisterPhysics( new A01EMPhysics("standard EM"));

  // Muon Physics ( Apply related processes to mu and tau
  RegisterPhysics(  new A01MuonPhysics("muon"));

   // Hadron Physics ( Apply related processes to hadrons )
  RegisterPhysics(  new A01HadronPhysics("hadron"));
// We do not use hadronic lists since v7.
  //RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron"));
  //RegisterPhysics(  new HadronPhysicsQGSP_BIC("hadron"));

  // Ion Physics ( Apply related processes to ions )
  RegisterPhysics( new A01IonPhysics("ion"));

}

A01PhysicsList::~A01PhysicsList()
{
}

void A01PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}



