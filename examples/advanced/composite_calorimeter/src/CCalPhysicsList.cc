///////////////////////////////////////////////////////////////////////////////
// File: CCalPhysicsList.cc
// Description: CCalPhysicsList provides Physics list for the simulation
///////////////////////////////////////////////////////////////////////////////

#include "CCalPhysicsList.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

#include "CCalProdGeneralPhysics.hh"
#include "CCalProdEMPhysics.hh"
#include "CCalProdMuonPhysics.hh"
#include "CCalProdHadronPhysics.hh"
#include "CCalProdIonPhysics.hh"

CCalPhysicsList::CCalPhysicsList(): G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics (new CCalProdGeneralPhysics("general") );

  // EM Physics
  RegisterPhysics (new CCalProdEMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics (new CCalProdMuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics (new CCalProdHadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics (new CCalProdIonPhysics("ion"));
}

CCalPhysicsList::~CCalPhysicsList() {}


void CCalPhysicsList::SetCuts() {
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}


