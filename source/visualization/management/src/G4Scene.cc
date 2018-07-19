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
// $Id: G4Scene.cc 101958 2016-12-12 08:04:35Z gcosmo $
//
// 
// Scene data  John Allison  19th July 1996.

#include "G4Scene.hh"

#include "G4Vector3D.hh"
#include "G4BoundingSphereScene.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4TransportationManager.hh"

#include <set>

G4Scene::G4Scene (const G4String& name):
  fName (name),
  fRefreshAtEndOfEvent(true),
  fRefreshAtEndOfRun(true),
  fMaxNumberOfKeptEvents(100)
{} // Note all other data members have default initial values.

G4Scene::~G4Scene () {}

G4bool G4Scene::AddRunDurationModel (G4VModel* pModel, G4bool warn)
{
  std::vector<Model>::const_iterator i;
  for (i = fRunDurationModelList.begin ();
       i != fRunDurationModelList.end (); ++i) {
    if (pModel -> GetGlobalDescription () ==
	i->fpModel->GetGlobalDescription ()) break;
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

  for (i = fRunDurationModelList.begin ();
       i != fRunDurationModelList.end (); ++i) {
    if (pModel -> GetGlobalTag () ==
        i->fpModel->GetGlobalTag ()) break;
  }
  if (i != fRunDurationModelList.end ()) {
    if (warn) {
      G4cout
      << "G4Scene::AddRunDurationModel: The tag \""
      << pModel->GetGlobalTag()
      << "\"\n  duplicates one already in scene \""
      << fName
      <<
  "\".\n  This may be intended but if not, you may inspect the scene with"
  "\n  \"/vis/scene/list\" and deactivate unwanted models with"
  "\n  \"/vis/scene/activateModel\". Or create a new scene."
      << G4endl;
    }
  }

  fRunDurationModelList.push_back (Model(pModel));

  CalculateExtent ();

  return true;
}

void G4Scene::CalculateExtent ()
{
  G4BoundingSphereScene boundingSphereScene;

  for (size_t i = 0; i < fRunDurationModelList.size(); i++) {
    if (fRunDurationModelList[i].fActive) {
      G4VModel* model = fRunDurationModelList[i].fpModel;
      if (model -> Validate()) {  // Validates and also recomputes extent.
	const G4VisExtent& thisExtent = model -> GetExtent ();
	G4double thisRadius = thisExtent.GetExtentRadius ();
	if (thisRadius > 0.) {
	  G4Point3D thisCentre = thisExtent.GetExtentCentre ();
	  thisCentre.transform (model -> GetTransformation ());
	  boundingSphereScene.AccrueBoundingSphere (thisCentre, thisRadius);
	}
      } else {
	G4ExceptionDescription ed;
	ed << "Invalid model \"" << model->GetGlobalDescription()
	   << "\".\n  Not included in extent calculation.";
	G4Exception
	  ("G4Scene::CalculateExtent",
	   "visman0201", JustWarning, ed);
      }
    }
  }

  for (size_t i = 0; i < fEndOfEventModelList.size(); i++) {
    if (fEndOfEventModelList[i].fActive) {
      G4VModel* model = fEndOfEventModelList[i].fpModel;
      if (model -> Validate()) {  // Validates and also recomputes extent.
	const G4VisExtent& thisExtent = model -> GetExtent ();
	G4double thisRadius = thisExtent.GetExtentRadius ();
	if (thisRadius > 0.) {
	  G4Point3D thisCentre = thisExtent.GetExtentCentre ();
	  thisCentre.transform (model -> GetTransformation ());
	  boundingSphereScene.AccrueBoundingSphere (thisCentre, thisRadius);
	}
      } else {
	G4ExceptionDescription ed;
	ed << "Invalid model \"" << model->GetGlobalDescription()
	   << "\".\n  Not included in extent calculation.";
	G4Exception
	  ("G4Scene::CalculateExtent",
	   "visman0201", JustWarning, ed);
      }
    }
  }

  for (size_t i = 0; i < fEndOfRunModelList.size(); i++) {
    if (fEndOfRunModelList[i].fActive) {
      G4VModel* model = fEndOfRunModelList[i].fpModel;
      if (model -> Validate()) {  // Validates and also recomputes extent.
	const G4VisExtent& thisExtent = model -> GetExtent ();
	G4double thisRadius = thisExtent.GetExtentRadius ();
	if (thisRadius > 0.) {
	  G4Point3D thisCentre = thisExtent.GetExtentCentre ();
	  thisCentre.transform (model -> GetTransformation ());
	  boundingSphereScene.AccrueBoundingSphere (thisCentre, thisRadius);
	}
      } else {
	G4ExceptionDescription ed;
	ed << "Invalid model \"" << model->GetGlobalDescription()
	   << "\".\n  Not included in extent calculation.";
	G4Exception
	  ("G4Scene::CalculateExtent",
	   "visman0201", JustWarning, ed);
      }
    }
  }

  fExtent = boundingSphereScene.GetBoundingSphereExtent ();
  fStandardTargetPoint = fExtent.GetExtentCentre ();
  if (fExtent.GetExtentRadius() <= 0.) {
	G4Exception
	  ("G4Scene::CalculateExtent",
	   "visman0202", JustWarning,
	   "Scene has no extent.  Please activate or add something."
	   "\nThe camera needs to have something to point at!"
           "\nAdd a volume. (You may need \"/run/initialize\".)"
           "\nOr use \"/vis/scene/add/extent\"."
	   "\n\"/vis/scene/list\" to see list of models.");
  }
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
		" SetVisAttributes (G4VisAttributes::GetInvisible());"
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
	fEndOfEventModelList[i].fpModel -> GetGlobalDescription ()) break;
  }
  if (i < nModels) {
    if (warn) {
      G4cout << "G4Scene::AddEndOfEventModel: a model \""
	     << pModel -> GetGlobalDescription ()
	     << "\"\n  is already in the end-of-event list of scene \""
	     << fName << "\"."
	     << G4endl;
    }
    return false;
  }
  fEndOfEventModelList.push_back (Model(pModel));
  return true;
}

G4bool G4Scene::AddEndOfRunModel (G4VModel* pModel, G4bool warn) {
  G4int i, nModels = fEndOfRunModelList.size ();
  for (i = 0; i < nModels; i++) {
    if (pModel -> GetGlobalDescription () ==
	fEndOfRunModelList[i].fpModel -> GetGlobalDescription ()) break;
  }
  if (i < nModels) {
    if (warn) {
      G4cout << "G4Scene::AddEndOfRunModel: a model \""
	     << pModel -> GetGlobalDescription ()
	     << "\"\n  is already in the end-of-run list of scene \""
	     << fName << "\"."
	     << G4endl;
    }
    return false;
  }
  fEndOfRunModelList.push_back (pModel);
  return true;
}

std::ostream& operator << (std::ostream& os, const G4Scene& scene) {

  size_t i;

  os << "Scene data:";

  os << "\n  Run-duration model list:";
  for (i = 0; i < scene.fRunDurationModelList.size (); i++) {
    if (scene.fRunDurationModelList[i].fActive) os << "\n  Active:   ";
    else os << "\n  Inactive: ";
    os << *(scene.fRunDurationModelList[i].fpModel);
  }

  os << "\n  End-of-event model list:";
  for (i = 0; i < scene.fEndOfEventModelList.size (); i++) {
    if (scene.fEndOfEventModelList[i].fActive) os << "\n  Active:   ";
    else os << "\n  Inactive: ";
    os << *(scene.fEndOfEventModelList[i].fpModel);
  }

  os << "\n  End-of-run model list:";
  for (i = 0; i < scene.fEndOfRunModelList.size (); i++) {
    if (scene.fEndOfRunModelList[i].fActive) os << "\n  Active:   ";
    else os << "\n  Inactive: ";
    os << *(scene.fEndOfRunModelList[i].fpModel);
  }

  os << "\n  Extent or bounding box: " << scene.fExtent;

  os << "\n  Standard target point:  " << scene.fStandardTargetPoint;

  os << "\n  End of event action set to \"";
  if (scene.fRefreshAtEndOfEvent) os << "refresh\"";
  else {
    os << "accumulate (maximum number of kept events: ";
    if (scene.fMaxNumberOfKeptEvents >= 0) os << scene.fMaxNumberOfKeptEvents;
    else os << "unlimited";
    os << ")";
  }

  os << "\n  End of run action set to \"";
  if (scene.fRefreshAtEndOfRun) os << "refresh";
  else os << "accumulate";
  os << "\"";

  return os;
}

G4bool G4Scene::operator != (const G4Scene& scene) const {
  if (
      (fRunDurationModelList.size () !=
       scene.fRunDurationModelList.size ())                 ||
      (fEndOfEventModelList.size () !=
       scene.fEndOfEventModelList.size ())                  ||
      (fEndOfRunModelList.size () !=
       scene.fEndOfRunModelList.size ())                    ||
      (fExtent               != scene.fExtent)              ||
      !(fStandardTargetPoint == scene.fStandardTargetPoint) ||
      fRefreshAtEndOfEvent   != scene.fRefreshAtEndOfEvent  ||
      fRefreshAtEndOfRun     != scene.fRefreshAtEndOfRun    ||
      fMaxNumberOfKeptEvents != scene.fMaxNumberOfKeptEvents
      ) return true;

  /* A complete comparison should, perhaps, include a comparison of
     individual models, but it is not easy to implement operator!= for
     all models.  Also, it would be unfeasible to ask users to
     implement opeerator!= if we ever get round to allowing
     user-defined models.  Moreover, there is no editing of G4Scene
     objects, apart from changing fRefreshAtEndOfEvent, etc; as far as
     models are concerned, all you can ever do is add them, so a test
     on size (above) is enough.

  for (size_t i = 0; i < fRunDurationModelList.size (); i++) {
    if (fRunDurationModelList[i] != scene.fRunDurationModelList[i])
      return true;
  }

  for (size_t i = 0; i < fEndOfEventModelList.size (); i++) {
    if (fEndOfEventModelList[i] != scene.fEndOfEventModelList[i])
      return true;
  }

  for (size_t i = 0; i < fEndOfRunModelList.size (); i++) {
    if (fEndOfRunModelList[i] != scene.fEndOfRunModelList[i])
      return true;
  }
  */

  return false;
}
