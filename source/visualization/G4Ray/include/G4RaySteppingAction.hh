// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RaySteppingAction.hh,v 2.0 1998/07/02 16:42:50 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// 16/Apr/1997  J. Allison:  For visualization/test/test19.
// moved here and renamed 14/April/1997 

#ifndef MYSTEPPINGACTION_HH
#define MYSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4Normal3D.hh"

class G4RaySteppingAction: public G4UserSteppingAction {
public:
  G4RaySteppingAction();
  void UserSteppingAction();
  G4bool WeHaveContact () const {return weHaveContact;}
  G4Normal3D GetSurfaceNormal () const {return normal;}
  const G4VisAttributes* GetVisAttributes () const {return pVisAttribs;}
  G4VPhysicalVolume* GetPhysicalVolume () const {return pVPV;}
  G4LogicalVolume* GetLogicalVolume () const {return pLV;}
  G4VSolid* GetSolid () const {return pSol;}
private:
  G4bool weHaveContact;
  G4Normal3D normal;
  const G4VisAttributes* pVisAttribs;
  G4VPhysicalVolume* pVPV;
  G4LogicalVolume* pLV;
  G4VSolid* pSol;
};

#endif
