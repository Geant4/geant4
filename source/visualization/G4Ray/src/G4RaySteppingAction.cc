// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RaySteppingAction.cc,v 2.1 1998/07/13 17:11:12 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// 16/Apr/1997  J. Allison:  For visualization/test/test19.

#include "G4RaySteppingAction.hh"

#include "G4SteppingManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"

G4RaySteppingAction::G4RaySteppingAction() :  weHaveContact(false), pVisAttribs(0){}

void G4RaySteppingAction::UserSteppingAction() {

  G4Step* piStep = GetSteppingManager () -> GetStep ();
  // G4cout <<"G4Step* piStep: "<< (void *) piStep << endl;
 
  G4Track* pTrack = piStep -> GetTrack();
  //G4cout << "G4Track* pTrack: "<<(void *) pTrack << endl;

  pVPV = pTrack -> GetNextVolume();
  //G4cout << "G4VPhysicalVolume* pVPV: "<<(void *) pVPV << endl;
       
  if (pVPV) {

    weHaveContact = true;

    pLV = pVPV -> GetLogicalVolume();
    //G4cout <<"G4LogicalVolume* pLV: "<< (void *) pLV << endl;
 
    pSol = pLV -> GetSolid();
    //G4cout <<"G4VSolid* pSol: "<<  (void *) pSol << endl;
   
    G4Point3D point = piStep -> GetPostStepPoint () -> GetPosition ();
    //G4cout << "Point: " << point << endl;
  
    G4Vector3D translation = pVPV -> GetTranslation();
    G4RotationMatrix* rotation = pVPV -> GetRotation();
    if (rotation) {
       point.transform(G4Transform3D(*rotation, -translation));
	        }
    else {
       point.transform(G4Translate3D(-translation));
         }
    //G4cout << "Point after transformation: " << point << endl;
 
    normal = pSol -> SurfaceNormal(point);
    //G4cout << "normal: " << normal << endl;

    pVisAttribs = pLV -> GetVisAttributes ();
    //G4cout << "Vis attributes: " << *pVisAttribs << endl;
  }
  else {
    weHaveContact = false;
  }

  pTrack -> SetTrackStatus(fStopAndKill);

}




































