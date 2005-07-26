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

#include "G4HumanPhantomSteppingAction.hh"
#include "G4HumanPhantomEventAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#ifdef G4ANALYSIS_USE
#include "G4HumanPhantomAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HumanPhantomSteppingAction::G4HumanPhantomSteppingAction(G4HumanPhantomEventAction* eventAction):event(eventAction)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HumanPhantomSteppingAction::~G4HumanPhantomSteppingAction()
{}


void G4HumanPhantomSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

  G4String name = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  //  Give Segmentation Fault... why?
  // G4String name2 =   aStep->GetTrack()->GetNextVolume()->GetName();

  //G4cout << name2 << G4endl;

  //if ((aStep->GetTrack()->GetNextVolume()->GetName() != "OutOfThisWorld"))
  //{
     G4double stepLength = aStep -> GetTrack() -> GetStepLength(); 
     event -> SetPath(stepLength);

     G4ThreeVector particlePosition = aStep -> GetTrack() -> GetPosition();

#ifdef G4ANALYSIS_USE

     if (particlePosition[2]/cm >= 0.) // section at the z-coord
      { 
	G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();  	
	analysis -> particleProjectionXY(particlePosition[0]/cm,particlePosition[1]/cm);
      }

     if (particlePosition[0]/cm >= 0.) // section at the x-coord
      { 
	G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
	analysis -> particleProjectionYZ(particlePosition[1]/cm,particlePosition[2]/cm);
      }

     if (particlePosition[1]/cm >= 50. and particlePosition[1]/cm <= 51.) // section at the y-coord
      { 
	G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
	analysis -> particleProjectionZX(particlePosition[2]/cm,particlePosition[0]/cm);
      }

#endif  

     // }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

