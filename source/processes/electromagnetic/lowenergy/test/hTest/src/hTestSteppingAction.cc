// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestSteppingAction -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestSteppingAction.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestSteppingAction::hTestSteppingAction(hTestRunAction* ra,
                                         hTestEventAction* ea,
                                         hTestDetectorConstruction* det)
 :detector(det),
  eventaction(ea),
  runaction(ra),
  trIDnow(-1),
  trIDold(-1),
  evnOld(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestSteppingAction::~hTestSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

  G4int evno = eventaction->GetEventNo(); 
  trIDnow = aStep->GetTrack()->GetTrackID();
  G4double tkin = aStep->GetTrack()->GetKineticEnergy(); 
  G4double theta = (aStep->GetTrack()->GetMomentumDirection()).theta();

  G4bool stop = false;
  G4bool primary = false;
  G4bool outAbs = false;

  if(tkin == 0.0) stop = true;
  if(0 == aStep->GetTrack()->GetParentID()) primary = true;
  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber" &&
     aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="World")
     outAbs = true;

  // new particle
  if(trIDnow != trIDold || evno != evnOld) {
    trIDold = trIDnow;
    evnOld = evno;
  }

  // After step in absorber
  if(outAbs) {
    if(theta < M_PI * 0.5) eventaction->AddLeakEnergy(tkin);
    else                   eventaction->AddBackEnergy(tkin);

    // primary particles
    if(primary) {
      runaction->SaveToTuple("TEND",tkin/MeV);
      runaction->SaveToTuple("TETA",theta/deg);      
    }
  }

  // Primary particle stop 
  if(primary && stop) {

    G4double xend = aStep->GetPostStepPoint()->GetPosition().x();
    G4double yend = aStep->GetPostStepPoint()->GetPosition().y();
    G4double zend = aStep->GetPostStepPoint()->GetPosition().z();
    runaction->SaveToTuple("XEND",xend/mm);      
    runaction->SaveToTuple("YEND",yend/mm);      
    runaction->SaveToTuple("ZEND",zend/mm);      
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

