#include "LXePhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "LXeGeneralPhysics.hh"
#include "LXeEMPhysics.hh"
#include "LXeMuonPhysics.hh"
#include "LXeOpticalPhysics.hh"

LXePhysicsList::LXePhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  // SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new LXeGeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new LXeEMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new LXeMuonPhysics("muon"));

   // Optical Physics
  RegisterPhysics(  new LXeOpticalPhysics("optical"));

}

LXePhysicsList::~LXePhysicsList()
{
}

void LXePhysicsList::SetCuts(){
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}




