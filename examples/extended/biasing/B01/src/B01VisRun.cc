#include "B01VisRun.hh"
#include "G4RunManager.hh"
#include "B01PhysicsList.hh"
#include "B01PrimaryGeneratorAction.hh"
#include "B01DetectorConstruction.hh"
#include "B01VisEventAction.hh"
#include "B01RunAction.hh"



//#include "G4DosimPhysics.hh"
//#include "B01SteppingAction.hh"

B01VisRun::B01VisRun()
  :
  fRunManager(new G4RunManager)
{
  if (!fRunManager) {
    G4std::G4Exception("B01VisRun::B01VisRun: new failed to create G4RunManager!");
  }
  fRunManager->SetUserAction(new B01VisEventAction);
  fRunManager->SetUserAction(new B01RunAction);
}
B01VisRun::~B01VisRun(){}

void B01VisRun::SetDetector(G4VPhysicalVolume &worldvol){

  
  fRunManager->
    SetUserInitialization(new B01DetectorConstruction(worldvol));
  
  fRunManager->SetUserInitialization(new B01PhysicsList);
  /*
  PhysicsList *ph = new PhysicsList;
  fRunManager->SetUserInitialization(ph);
  ph->ChangeModel(LEinelastic);
  fRunManager->SetUserAction(new B01SteppingAction);
  */


  fRunManager->SetUserAction(new B01PrimaryGeneratorAction);  

}
void B01VisRun::Initialize(){
  fRunManager->Initialize();
}


