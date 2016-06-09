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
// $Id: GammaRayTelPhysicsList.cc,v 1.5 2003/06/16 16:46:30 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
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

#include "GammaRayTelGeneralPhysics.hh"
#include "GammaRayTelEMPhysics.hh"
#include "GammaRayTelMuonPhysics.hh"
#include "GammaRayTelHadronPhysics.hh"
#include "GammaRayTelIonPhysics.hh"

GammaRayTelPhysicsList::GammaRayTelPhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new GammaRayTelGeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new GammaRayTelEMPhysics("electromagnetic"));

  // Muon Physics
  RegisterPhysics(  new GammaRayTelMuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new GammaRayTelHadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new GammaRayTelIonPhysics("ion"));


}

GammaRayTelPhysicsList::~GammaRayTelPhysicsList()
{
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
  SetCutsWithDefault();

 
  // Production thresholds for detector regions

  G4String regName[] = {"Calorimeter","Tracker"};
  G4double cutValue[] = {1*mm, 0.1*mm};
  
  for(G4int i=0;i<2;i++)
  { 
    G4Region* reg = G4RegionStore::GetInstance()->GetRegion(regName[i]);
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(cutValue[i]);
    reg->SetProductionCuts(cuts);
  }
}




