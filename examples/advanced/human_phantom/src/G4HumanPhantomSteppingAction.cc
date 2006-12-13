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
