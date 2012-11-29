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
/// \file runAndEvent/RE02/src/RE02PhysicsList.cc
/// \brief Implementation of the RE02PhysicsList class
//
// $Id$
// --------------------------------------------------------------
//
// 28-Jan-04 Add QGSP_BERT and QGSP_BIC for hadronic lists. T. Koi
// 22-Nov-04 Comment out QGSP_BERT and QGSP_BIC
//           Output Notificaiton message             
//           All Particles are created in GeneralPhysics 

#include "RE02PhysicsList.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
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

#include "RE02GeneralPhysics.hh"
#include "RE02EMPhysics.hh"
#include "RE02MuonPhysics.hh"
#include "RE02HadronPhysics.hh"
#include "RE02IonPhysics.hh"

//#include "HadronPhysicsQGSP_BERT.hh"
//#include "HadronPhysicsQGSP_BIC.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02PhysicsList::RE02PhysicsList() : G4VModularPhysicsList()
{

  G4cout
    << "You are using the RE02PhysicsList\n"
    << "This PhysicsList originally comes from example/extended/analysis/A01,\n"
    << "and is modified in Hadron Physics in order to involve Binary Cascade\n"
    << "at low energy region and inelastic process for generic ions.\n"
    << "Full set of particles (barions bosons and mesons) will be created and\n"
    << "Standard EM Physics and Low & High Energy parameterized models will be "
    << "applied.\n"
    << "RE02PhysicsList is optimized for robustness and not for any particular "
    << "usage.\n"
    << "For the hadronic physics, educated guesses of physics list are "
    << "prepared\n"
    << "for various use cases.\n"
    << "When you will start REAL calculations for your own interest, please\n"
    << "consider the usage of hadronic_lists instead of RE02PhysicsLists.\n"
    << "More information can also be found from the Geant4 HyperNews.\n"
    << "http://geant4-hn.slac.stanford.edu:5090/Geant4-HyperNews/index\n"
    << G4endl;

  // default cut value  (1.0mm)
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics ( Create ALL Particle and apply Decay )
  RegisterPhysics( new RE02GeneralPhysics("general") );

  // EM Physics ( Apply related Processes to gamma and e-/+)
  RegisterPhysics( new RE02EMPhysics("standard EM"));

  // Muon Physics ( Apply related processes to mu and tau
  RegisterPhysics(  new RE02MuonPhysics("muon"));

   // Hadron Physics ( Apply related processes to hadrons )
   RegisterPhysics(  new RE02HadronPhysics("hadron"));
   // We do not use hadronic lists since v7.
  //RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron"));
  //RegisterPhysics(  new HadronPhysicsQGSP_BIC("hadron"));

  // Ion Physics ( Apply related processes to ions )
  RegisterPhysics( new RE02IonPhysics("ion"));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02PhysicsList::~RE02PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE02PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}



