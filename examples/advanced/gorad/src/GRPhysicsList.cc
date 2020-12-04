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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRPhysicsList.cc
//   Gorad Physics List
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRPhysicsList.hh"
#include "GRPhysicsListMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysListFactory.hh"

GRPhysicsList::GRPhysicsList()
: PLName("FTFP_BERT"), physList(nullptr),
  EM_opt("Op_0"), Had_opt("FTFP_BERT"),
  addHP(false), addRDM(false), addRMC(false), addOptical(false),
  stepLimit_opt(-1)
{
  factory = nullptr;
  messenger = new GRPhysicsListMessenger(this);

  globalCuts[0] = 0.7*mm; // e-
  globalCuts[1] = 0.7*mm; // e+
  globalCuts[2] = 0.7*mm; // gamma
  globalCuts[3] = 0.7*mm; // proton
}

GRPhysicsList::~GRPhysicsList()
{
  delete factory;
  if(physList) delete physList;
  delete messenger;
}

void GRPhysicsList::ConstructParticle()
{
  if(!physList) GeneratePL();
  physList->ConstructParticle();
}

void GRPhysicsList::ConstructProcess()
{
  if(!physList) GeneratePL();
  physList->ConstructProcess();
}

#include "G4Region.hh"
#include "G4ProductionCuts.hh"

void GRPhysicsList::SetCuts()
{
  if(!physList) GeneratePL();
  physList->SetCutValue(globalCuts[2],"gamma"); // gamma should be defined first!
  physList->SetCutValue(globalCuts[0],"e-");
  physList->SetCutValue(globalCuts[1],"e+");
  physList->SetCutValue(globalCuts[3],"proton");
}
  
void GRPhysicsList::SetGlobalCuts(G4double val)
{
  for(G4int i=0; i<4; i++)
  { SetGlobalCut(i,val); }
  if(physList) SetCuts();
}

void GRPhysicsList::SetGlobalCut(G4int i, G4double val) 
{
  globalCuts[i] = val; 
  if(physList) SetCuts();
}

void GRPhysicsList::GeneratePLName()
{
  G4String plname = Had_opt;
  if(addHP && Had_opt != "Shielding") plname += "_HP";

  G4String EMopt = "";
  if(EM_opt=="Op_1") EMopt = "_EMV"; 
  else if(EM_opt=="Op_3") EMopt = "_EMY"; 
  else if(EM_opt=="Op_4") EMopt = "_EMZ"; 
  else if(EM_opt=="LIV")  EMopt = "_LIV"; 
  else if(EM_opt=="LIV_Pol") G4cout << "EM option <LIV_Pol> is under development." << G4endl;
  plname += EMopt;

  auto valid = factory->IsReferencePhysList(plname);
  if(valid) 
  { PLName = plname; }
  else
  {
    G4ExceptionDescription ed;
    ed << "Physics List <" << plname << "> is not a valid reference physics list.";
    G4Exception("GRPhysicsList::GeneratePLName()","GRPHYS0001",
                FatalException,ed);
  }
}

#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4GenericBiasingPhysics.hh"
#include "GRParallelWorldPhysics.hh"
#include "G4ProcessTable.hh"
#include "G4EmParameters.hh"
#include "G4HadronicParameters.hh"

void GRPhysicsList::GeneratePL()
{
  if(physList) return;

  factory = new G4PhysListFactory();
  GeneratePLName();
  physList = factory->GetReferencePhysList(PLName);
  G4cout << "Creating " << PLName << " physics list ################ " << applyGeomImpBias << G4endl; 

  if(addRDM && Had_opt != "Shielding")
  { physList->RegisterPhysics(new G4RadioactiveDecayPhysics()); 
    G4cout << "Adding G4RadioactiveDecayPhysics ################ " << G4endl; }

  if(addRMC)
  { G4cout << "Reverse Monte Calro option is under development." << G4endl; }

  if(stepLimit_opt>=0)
  { physList->RegisterPhysics(new G4StepLimiterPhysics()); 
    G4cout << "Adding G4StepLimiterPhysics ################ " << G4endl; }

  if(addOptical) // Optical processes should come last!
  { physList->RegisterPhysics(new G4OpticalPhysics()); 
    G4cout << "Adding G4OpticalPhysics ################ " << G4endl; }

  if(applyGeomImpBias) // Geometry Importance Biasing with parallel world
  {
    physList->RegisterPhysics(new GRParallelWorldPhysics("GeomBias",false));
    G4cout << "Adding G4GenericBiasingPhysics for GeomBias ################ " << G4endl;
  }

  G4int verbose = G4ProcessTable::GetProcessTable()->GetVerboseLevel();
  physList->SetVerboseLevel(verbose);
  G4EmParameters::Instance()->SetVerbose(verbose);
  G4HadronicParameters::Instance()->SetVerboseLevel(verbose);
}

#include "G4RegionStore.hh"

G4Region* GRPhysicsList::FindRegion(const G4String& reg) const
{
  auto store = G4RegionStore::GetInstance();
  return store->GetRegion(reg);
}

G4Region* GRPhysicsList::SetLocalCut(const G4String& reg,G4int i,G4double val)
{
  auto regPtr = FindRegion(reg);
  if(!regPtr) return regPtr;

  auto cuts = regPtr->GetProductionCuts();
  if(!cuts) 
  {
    cuts = new G4ProductionCuts();
    regPtr->SetProductionCuts(cuts);
  }

  cuts->SetProductionCut(val,i);
  return regPtr;
}

G4double GRPhysicsList::GetLocalCut(const G4String& reg,G4int i) const
{
  auto regPtr = FindRegion(reg);
  G4double val = -1.0;
  if(regPtr) 
  { 
    auto cuts = regPtr->GetProductionCuts();
    if(cuts) val = cuts->GetProductionCut(i);
  }
  return val;
}

#include "G4UserLimits.hh"

G4Region* GRPhysicsList::SetLocalStepLimit(const G4String& reg,G4double val)
{
  auto regPtr = FindRegion(reg);
  if(!regPtr) return regPtr;

  auto uLim = regPtr->GetUserLimits();
  if(!uLim)
  {
    uLim = new G4UserLimits(val);
    regPtr->SetUserLimits(uLim);
  }
  else
  { uLim->SetMaxAllowedStep(val); }
  return regPtr;
}

#include "G4Track.hh"
G4double GRPhysicsList::GetLocalStepLimit(const G4String& reg) const
{
  static G4Track dummyTrack;
  auto regPtr = FindRegion(reg);
  G4double val = -1.0;
  if(regPtr)
  {
    auto uLim = regPtr->GetUserLimits();
    if(uLim) val = uLim->GetMaxAllowedStep(dummyTrack);
  }
  return val;
}

void GRPhysicsList::SetGlobalStepLimit(G4double val)
{ SetLocalStepLimit("DefaultRegionForTheWorld",val); }

G4double GRPhysicsList::GetGlobalStepLimit() const
{ return GetLocalStepLimit("DefaultRegionForTheWorld"); }


