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
// $Id: GammaRayTelPhysicsList.cc,v 1.8 2009-11-18 15:59:05 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "GammaRayTelPhysicsList.hh"

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

#include "GammaRayTelParticles.hh"
#include "GammaRayTelGeneralPhysics.hh"
#include "GammaRayTelEMstdPhysics.hh"
#include "GammaRayTelEMlowePhysics.hh"
#include "GammaRayTelMuonPhysics.hh"
#include "GammaRayTelHadronPhysics.hh"
#include "GammaRayTelIonPhysics.hh"

#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4EmProcessOptions.hh"


#include "GammaRayTelPhysicsListMessenger.hh"

GammaRayTelPhysicsList::GammaRayTelPhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  pMessenger = new GammaRayTelPhysicsListMessenger(this);

  // Particles

 
  RegisterPhysics( new GammaRayTelParticles("particles") );

  G4cout << "PARTICLES DONE" << G4endl;

  // EM physics

  emPhysicsList = new GammaRayTelEMstdPhysics;
  emName = G4String("Standard EM");


  // General Physics

  RegisterPhysics( new GammaRayTelGeneralPhysics("general") );

  G4cout << "GENERAL DONE" << G4endl;


  // Muon Physics

  RegisterPhysics(  new GammaRayTelMuonPhysics("muon"));


  G4cout << "MUON DONE" << G4endl;

   // Hadron Physics
  RegisterPhysics(  new GammaRayTelHadronPhysics("hadron"));

  G4cout << "HADRONS DONE" << G4endl;

  // Ion Physics
  RegisterPhysics( new GammaRayTelIonPhysics("ion"));


  G4cout << "IONS DONE" << G4endl;

}

GammaRayTelPhysicsList::~GammaRayTelPhysicsList()
{
  delete pMessenger;
  delete emPhysicsList;
}

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

void GammaRayTelPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "GammaRayTelPhysicsList::SetCuts: default cut length : "
         << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  

  // These values are used as the default production thresholds
  // for the world volume.

  G4cout << "CUT STD" << G4endl;

  SetCutsWithDefault();


}

void GammaRayTelPhysicsList::SetRegionCut(G4double cutvalue)
{

 SetCutsWithDefault();

  if (verboseLevel >0){
    G4cout << "GammaRayTelPhysicsList::SetCuts: default cut length : "
         << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  
  G4cout << "CUTS NEW" << G4endl;
  
  // Production thresholds for detector regions

  G4String regName[] = {"Calorimeter","Tracker"};
  //  G4double cutValue[] = {1*mm, 0.1*mm};
  G4double cutValue[] = {cutvalue, cutvalue/10.};


  
  for(G4int i=0;i<2;i++)
  { 
    G4Region* reg = G4RegionStore::GetInstance()->GetRegion(regName[i]);
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(cutValue[i]);
    reg->SetProductionCuts(cuts);
  }

}


/////////////////////////////////////////////////////////////////////////////

void GammaRayTelPhysicsList::AddPhysicsList(const G4String& name)
{

  if (verboseLevel>1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == "Standard EM") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new GammaRayTelEMstdPhysics();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: EM Standard" << G4endl;
    
  } else if (name == "LowE EM") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new GammaRayTelEMlowePhysics();;
     G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: EM LowE" << G4endl;
  } 
  else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }


  G4cout << "REGISTRATION DONE " << G4endl;

    RegisterPhysics(emPhysicsList);

}
