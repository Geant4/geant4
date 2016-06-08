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
// $Id: G4VisCommandsCamera.cc,v 1.7.4.1 2001/06/28 19:16:15 gunter Exp $
// GEANT4 tag $Name:  $
//
// 

#include "G4VisCommandsCamera.hh"

#include "G4UnitsTable.hh"

void G4VisCommandCameraReset::SetValue () {
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  if (pVMan -> IsValidView ()) {
    const G4Scene* pScene = pVMan -> GetCurrentScene ();
    G4ViewParameters& vp = pVMan -> SetCurrentViewParameters ();
    vp.SetCurrentTargetPoint (G4Point3D());
    vp.SetZoomFactor (1.);
    vp.SetDolly (0.);
    vp.SetViewpointDirection (G4Vector3D (0., 0., 1.));
    vp.SetUpVector (G4Vector3D (0., 1., 0.));
    G4cout << "Target point reset to centre of scene, ("
	   << G4BestUnit (pScene -> GetStandardTargetPoint (), "Length")
	   << ")";
    G4cout << "\nZoom factor reset to 1.";
    G4cout << "\nDolly distance reset to 0.";
    G4cout << "\nViewpoint direction reset to +z.";
    G4cout << "\nUp vector reset to +y.";
    G4cout << G4endl;
    if (pVMan -> GetVerboseLevel () > 1) {
      pVMan -> PrintCurrentView ();
    }
  }
}
