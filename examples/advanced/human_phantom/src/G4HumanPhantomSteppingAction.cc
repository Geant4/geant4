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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4HumanPhantomSteppingAction.hh"
#include "G4HumanPhantomEventAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#ifdef G4ANALYSIS_USE
#include "G4HumanPhantomAnalysisManager.hh"
#endif
          
G4HumanPhantomSteppingAction::G4HumanPhantomSteppingAction()
{ }

G4HumanPhantomSteppingAction::~G4HumanPhantomSteppingAction()
{}

void G4HumanPhantomSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  
  G4String name = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

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

  if ((particlePosition[1]/cm >= 50.) && (particlePosition[1]/cm <= 51.)) 
    // section at the y-coord
    { 
      G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
      analysis -> particleProjectionZX(particlePosition[2]/cm,particlePosition[0]/cm);
    }

#endif
  
}
