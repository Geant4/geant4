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
//
// $Id: G4Scene.cc,v 1.24 2009-11-04 12:49:16 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Scene data  John Allison  19th July 1996.

#include "G4Scene.hh"

#include "G4Vector3D.hh"
#include "G4BoundingSphereScene.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4TransportationManager.hh"

G4Scene::G4Scene (const G4String& name):
  fName (name),
  fRefreshAtEndOfEvent(true),
  fRefreshAtEndOfRun(true),
  fMaxNumberOfKeptEvents(0)
{} // Note all other data members have default initial values.

G4Scene::~G4Scene () {}

G4bool G4Scene::AddRunDurationModel (G4VModel* pModel, G4bool warn) {
  std::vector<G4VModel*>::const_iterator i;
  for (i = fRunDurationModelList.begin ();
       i != fRunDurationModelList.end (); ++i) {
    if (pModel -> GetGlobalDescription () ==
	(*i) -> GetGlobalDescription ()) break;
  }
  if (i != fRunDurationModelList.end ()) {
    if (warn) {
      G4cout << "G4Scene::AddRunDurationModel: model \""
	     << pModel -> GetGlobalDescription ()
	     << "\"\n  is already in the run-duration list of scene \""
	     << fName
	     << "\"."
	     << G4endl;
    }
    return false;
  }
  fRunDurationModelList.push_back (pModel);
  CalculateExtent ();
  return true;
}

void G4Scene::CalculateExtent () {
  G4int nModels = fRunDurationModelList.size ();
  G4BoundingSphereScene boundingSphereScene;
  for (G4int i = 0; i < nModels; i++) {
    const G4VisExtent& thisExtent =
      fRunDurationModelList[i] -> GetExtent ();
    G4Point3D thisCentre = thisExtent.GetExtentCentre ();
    G4double thisRadius = thisExtent.GetExtentRadius ();
    thisCentre.transform (fRunDurationModelList[i] -> GetTransformation ());
    boundingSphereScene.AccrueBoundingSphere (thisCentre, thisRadius);
  }
  fExtent = boundingSphereScene.GetBoundingSphereExtent ();
  fStandardTargetPoint = fExtent.GetExtentCentre ();
}

G4bool G4Scene::AddWorldIfEmpty (G4bool warn) {
  G4bool successful = true;
  if (IsEmpty ()) {
    successful = false;
    G4VPhysicalVolume* pWorld =
      G4TransportationManager::GetTransportationManager ()
      -> GetNavigatorForTracking () -> GetWorldVolume ();
    if (pWorld) {
      const G4VisAttributes* pVisAttribs =
	pWorld -> GetLogicalVolume () -> GetVisAttributes ();
      if (!pVisAttribs || pVisAttribs -> IsVisible ()) {
	if (warn) {
	  G4cout << 
	    "Your \"world\" has no vis attributes or is marked as visible."
	    "\n  For a better view of the contents, mark the world as"
	    " invisible, e.g.,"
	    "\n  myWorldLogicalVol ->"
	    " SetVisAttributes (G4VisAttributes::Invisible);"
		 << G4endl;
	}
      }
      successful = AddRunDurationModel (new G4PhysicalVolumeModel (pWorld));
      // Note: default depth and no modeling parameters.
      if (successful) {
	if (warn) {
	  G4cout <<
    "G4Scene::AddWorldIfEmpty: The scene was empty of run-duration models."
    "\n  \"world\" has been added.";
	  G4cout << G4endl;
	}
      }
    }
  }
  return successful;
}

G4bool G4Scene::AddEndOfEventModel (G4VModel* pModel, G4bool warn) {
  G4int i, nModels = fEndOfEventModelList.size ();
  for (i = 0; i < nModels; i++) {
    if (pModel -> GetGlobalDescription () ==
	fEndOfEventModelList [i] -> GetGlobalDescription ()) break;
  }
  if (i < nModels) {
    delete fEndOfEventModelList[i];
    fEndOfEventModelList[i] = pModel;
    if (warn) {
      G4cout << "G4Scene::AddEndOfEventModel: a model \""
	     << pModel -> GetGlobalDescription ()
	     << "\"\n  is already in the end-of-event list of scene \""
	     << fName <<
	"\".\n  The old model has been deleted; this new model replaces it."
	     << G4endl;
    }
    return true;  // Model replaced sucessfully.
  }
  fEndOfEventModelList.push_back (pModel);
  return true;
}

G4bool G4Scene::AddEndOfRunModel (G4VModel* pModel, G4bool warn) {
  G4int i, nModels = fEndOfRunModelList.size ();
  for (i = 0; i < nModels; i++) {
    if (pModel -> GetGlobalDescription () ==
	fEndOfRunModelList [i] -> GetGlobalDescription ()) break;
  }
  if (i < nModels) {
    delete fEndOfRunModelList[i];
    fEndOfRunModelList[i] = pModel;
    if (warn) {
      G4cout << "G4Scene::AddEndOfRunModel: a model \""
	     << pModel -> GetGlobalDescription ()
	     << "\"\n  is already in the end-of-run list of scene \""
	     << fName <<
	"\".\n  The old model has been deleted; this new model replaces it."
	     << G4endl;
    }
    return true;  // Model replaced sucessfully.
  }
  fEndOfRunModelList.push_back (pModel);
  return true;
}

std::ostream& operator << (std::ostream& os, const G4Scene& s) {

  size_t i;

  os << "Scene data:";

  os << "\n  Run-duration model list:";
  for (i = 0; i < s.fRunDurationModelList.size (); i++) {
    os << "\n  " << *(s.fRunDurationModelList[i]);
  }

  os << "\n  End-of-event model list:";
  for (i = 0; i < s.fEndOfEventModelList.size (); i++) {
    os << "\n  " << *(s.fEndOfEventModelList[i]);
  }

  os << "\n  End-of-run model list:";
  for (i = 0; i < s.fEndOfRunModelList.size (); i++) {
    os << "\n  " << *(s.fEndOfRunModelList[i]);
  }

  os << "\n  Extent or bounding box: " << s.fExtent;

  os << "\n  Standard target point:  " << s.fStandardTargetPoint;

  os << "\n  End of event action set to \"";
  if (s.fRefreshAtEndOfEvent) os << "refresh\"";
  else {
    os << "accumulate (maximum number of kept events: ";
    if (s.fMaxNumberOfKeptEvents >= 0) os << s.fMaxNumberOfKeptEvents;
    else os << "unlimited";
    os << ")";
  }

  os << "\n  End of run action set to \"";
  if (s.fRefreshAtEndOfRun) os << "refresh";
  else os << "accumulate";
  os << "\"";

  return os;
}

G4bool G4Scene::operator != (const G4Scene& s) const {
  if (
      (fRunDurationModelList.size () !=
       s.fRunDurationModelList.size ())                 ||
      (fExtent               != s.fExtent)              ||
      !(fStandardTargetPoint == s.fStandardTargetPoint) ||
      fRefreshAtEndOfEvent   != s.fRefreshAtEndOfEvent  ||
      fRefreshAtEndOfRun     != s.fRefreshAtEndOfRun    ||
      fMaxNumberOfKeptEvents != s.fMaxNumberOfKeptEvents
      ) return true;

  for (size_t i = 0; i < fRunDurationModelList.size (); i++) {
    if (fRunDurationModelList[i] != s.fRunDurationModelList[i])
      return true;
  }

  return false;
}
