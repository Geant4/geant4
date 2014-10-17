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
// $Id$
// ------------------------------------------------------------
// Geant4 class implementation file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "AXPETSteppingAction.hh"

#include "AXPETDetectorConstruction.hh"
#include "AXPETRunAction.hh"

#include "G4Track.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ios.hh"
#include <iomanip>
#include "G4UImanager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

AXPETSteppingAction::AXPETSteppingAction(AXPETDetectorConstruction* DET)
  : detector (DET)
{
}


AXPETSteppingAction::~AXPETSteppingAction()
{;}


void AXPETSteppingAction::UserSteppingAction(const G4Step* aStep)
{

  // track informations
  //
  // const G4StepPoint* prePoint = aStep->GetPreStepPoint();   
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition(); 

  G4String partName  = particle->GetParticleName();
  G4VPhysicalVolume* postVolume_phys=endPoint->GetPhysicalVolume();
  G4String solidName=" out of Crystal ";
  if(postVolume_phys)    // Reach to out of the world
  solidName = postVolume_phys->GetLogicalVolume()->GetSolid()->GetName();
  //Detect optical photons going outside Volume

  if(partName=="opticalphoton")
  {
    // Additional Printing on each Step
    //
    // G4cout<<aStep->GetTrack()->GetCurrentStepNumber()<<" StartPoint="<<prePoint->GetPosition()
    // <<" EndPoint="<<endPoint->GetPosition()<<" NextSolid= "<<solidName<<G4endl;
    // G4cout<<aStep->GetTrack()->GetCurrentStepNumber()<<" StartDirection="<<prePoint->GetMomentumDirection()
    // <<" EndDirection="<<endPoint->GetMomentumDirection()<<" NextSolid= "<<solidName<<G4endl;
    G4double x=endPoint->GetPosition().x()/mm;
    G4double y=endPoint->GetPosition().y()/mm;
    //
    // Condition to detect Escaping Optical Photon: 
    //
    if(std::sqrt(x*x+y*y)>(detector->GetExtend()*2.)){
      G4cerr.precision(16);
      G4cerr << "ERROR - UserSteppingAction::OpticalPhoton is out of Crystal Solid" << G4endl
	     <<" Wrong Point is "<<xStep<<"  "<<yStep<<"  "<<zStep << G4endl
             <<" Wrong Direction is "<<xDirection<<"  "<<yDirection<<"  "<<zDirection<< G4endl;
      G4VSolid *solid=detector->GetSolid();
      G4ThreeVector point = G4ThreeVector(xStep,yStep,zStep);
      G4ThreeVector dir = G4ThreeVector(xDirection,yDirection,zDirection); 
        G4ThreeVector norm,*pNorm;
        G4bool *pgoodNorm, goodNorm, calcNorm=true;
        pNorm=&norm;
        pgoodNorm=&goodNorm;

       EInside surface=solid->Inside(point);
       if(surface == kInside){
	 G4cout<<"WrongPoint is Inside DistanceToOut="<<solid->DistanceToOut(point,dir,calcNorm,pgoodNorm,pNorm)<<" Norm(DistToOut)="<<norm<<G4endl;
       }
       else if(surface == kOutside){
       G4cout<<"WrongPoint is Outside DistanceToIn(p,+v)="<<solid->DistanceToIn(point,dir)<<G4endl;
       G4cout<<"                      DistanceToIn(p,-v)="<<solid->DistanceToIn(point,-dir)<<G4endl;
       }
       else{
       G4cout<<"WrongPoint is On Surface DistanceToIn(p,+v)="<<solid->DistanceToIn(point,dir)<<G4endl;
       G4cout<<"                         DistanceToIn(p,-v)="<<solid->DistanceToIn(point,-dir)<<G4endl;
       G4cout<<"                         DistanceToOut="<<solid->DistanceToOut(point,dir,calcNorm,pgoodNorm,pNorm)<<" Norm(DistanceToOut)="<<norm<<G4endl;
       G4cout<<"                         SurfaceNormal=="<<solid->SurfaceNormal(point)<<G4endl;
       }
       if(detector->GetAbortAction()){
       G4Exception( "AXPETSteppingAction", "AXPET001",
                  FatalException, "Optical Photon outside Crystal" );
       }
       else{
       G4Exception( "AXPETSteppingAction", "AXPET002",
                  JustWarning, "Optical Photon outside Crystal" );
       }
    }
    // Save values of Step
   
    xStep=x;
    yStep=y;
    zStep=endPoint->GetPosition().z()/mm;
    xDirection=endPoint->GetMomentumDirection().x();
    yDirection=endPoint->GetMomentumDirection().y();
    zDirection=endPoint->GetMomentumDirection().z();
  }
}
