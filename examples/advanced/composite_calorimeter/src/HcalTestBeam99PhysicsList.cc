///////////////////////////////////////////////////////////////////////////////
// File: HcalTestBeam99PhysicsList.cc
// Date: 04/2002 S.Banerjee
// Modifications: 
///////////////////////////////////////////////////////////////////////////////

#include "HcalTestBeam99PhysicsList.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

#include "CMSProdGeneralPhysics.hh"
#include "CMSProdEMPhysics.hh"
#include "CMSProdMuonPhysics.hh"
#include "CMSProdHadronPhysics.hh"
#include "CMSProdIonPhysics.hh"

HcalTestBeam99PhysicsList::HcalTestBeam99PhysicsList(): G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics (new CMSProdGeneralPhysics("general") );

  // EM Physics
  RegisterPhysics (new CMSProdEMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics (new CMSProdMuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics (new CMSProdHadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics (new CMSProdIonPhysics("ion"));
}

HcalTestBeam99PhysicsList::~HcalTestBeam99PhysicsList() {}


void HcalTestBeam99PhysicsList::SetCuts() {
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}


