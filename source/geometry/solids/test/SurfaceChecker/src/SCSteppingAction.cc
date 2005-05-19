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
// $Id: SCSteppingAction.cc,v 1.1 2005-05-19 13:07:29 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SCSteppingAction.hh"
#include "G4SteppingManager.hh"

#include "SCSurfacePoint.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCSteppingAction::SCSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

  G4Track* aTrack = aStep->GetTrack();

  G4StepPoint*  PreStepPoint= aStep->GetPreStepPoint(); 
  G4StepPoint*  PostStepPoint= aStep->GetPostStepPoint(); 
  G4ThreeVector PostHitPoint = PostStepPoint->GetPosition();
  G4ThreeVector PreHitPoint = PreStepPoint->GetPosition();

  G4cout.precision(16) ;

  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector m = (PostHitPoint - PreHitPoint).unit() ;

  G4ThreeVector truepoint = spoint.GetSurfacePoint() ;
  G4double delta = (truepoint - PostHitPoint).mag() ;

  G4cout << "Intersection "  
	 << PostHitPoint << " "  
	 << delta  <<  G4endl ;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

