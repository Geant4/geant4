// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Scene.cc,v 1.6 2001-02-23 15:43:22 johna Exp $
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
  fName (name)
{} // Note all data members have default initial values.

G4Scene::~G4Scene () {}

G4bool G4Scene::AddRunDurationModel (G4VModel* pModel) {
  G4int i, nModels = fRunDurationModelList.size ();
  for (i = 0; i < nModels; i++) {
    if (pModel -> GetGlobalDescription () ==
	fRunDurationModelList [i] -> GetGlobalDescription ()) break;
  }
  if (i < nModels) {
    G4cout << "G4Scene::AddRunDurationModel: model \""
	   << pModel -> GetGlobalDescription ()
	   << "\"\n  is already in the run-duration list of scene \""
	   << fName
	   << "\"."
	   << G4endl;
    return false;
  }
  fRunDurationModelList.push_back (pModel);
  nModels = fRunDurationModelList.size ();  // ...has increased by 1...
  G4BoundingSphereScene boundingSphereScene;
  for (i = 0; i < nModels; i++) {
    const G4VisExtent& thisExtent =
      fRunDurationModelList[i] -> GetExtent ();
    G4Point3D thisCentre = thisExtent.GetExtentCentre ();
    G4double thisRadius = thisExtent.GetExtentRadius ();
    thisCentre.transform (fRunDurationModelList[i] -> GetTransformation ());
    boundingSphereScene.AccrueBoundingSphere (thisCentre, thisRadius);
  }
  fExtent = boundingSphereScene.GetBoundingSphereExtent ();
  fStandardTargetPoint = fExtent.GetExtentCentre ();
  return true;
}

G4bool G4Scene::AddWorldIfEmpty () {
  G4bool successful = false;
  if (IsEmpty ()) {
    G4VPhysicalVolume* pWorld =
      G4TransportationManager::GetTransportationManager ()
      -> GetNavigatorForTracking () -> GetWorldVolume ();
    if (pWorld) {
      const G4VisAttributes* pVisAttribs =
	pWorld -> GetLogicalVolume () -> GetVisAttributes ();
      if (!pVisAttribs || pVisAttribs -> IsVisible ()) {
	G4cout << 
	  "Your \"world\" has no vis attributes or is marked as visible."
	  "\n  For a better view of the contents, mark the world as"
	  " invisible, e.g.,"
	  "\n  myWorldLogicalVol ->"
	  " SetVisAttributes (G4VisAttributes::Invisible);"
	       << G4endl;
      }
      successful = AddRunDurationModel (new G4PhysicalVolumeModel (pWorld));
      // Note: default depth and no modeling parameters.
      if (successful) {
	G4cout <<
	  "G4Scene::AddWorldIfEmpty: The scene was empty,"
	  "\n   \"world\" has been added.";
	G4cout << G4endl;
      }
    }
  }
  return successful;
}

G4bool G4Scene::AddEndOfEventModel (G4VModel* pModel) {
  G4int i, nModels = fEndOfEventModelList.size ();
  for (i = 0; i < nModels; i++) {
    if (pModel -> GetGlobalDescription () ==
	fEndOfEventModelList [i] -> GetGlobalDescription ()) break;
  }
  if (i < nModels) {
    G4cout << "G4Scene::AddEndOfEventModel: model \""
	   << pModel -> GetGlobalDescription ()
	   << "\"\n  is already in the run-duration list of scene \""
	   << fName
	   << "\"."
	   << G4endl;
    return false;
  }
  fEndOfEventModelList.push_back (pModel);
  return true;
}

void G4Scene::Clear () {
  int i;
  for (i = 0; i < fRunDurationModelList.size(); ++i) {
    delete fRunDurationModelList[i];
  }
  for (i = 0; i < fEndOfEventModelList.size(); ++i) {
    delete fEndOfEventModelList[i];
  }
}

G4std::ostream& operator << (G4std::ostream& os, const G4Scene& d) {

  os << "Scene data:";

  os << "\n  Run-duration model list:";
  for (int i = 0; i < d.fRunDurationModelList.size (); i++) {
    os << "\n  " << *(d.fRunDurationModelList[i]);
  }

  os << "\n  End-of-event model list:";
  for (int ii = 0; ii < d.fEndOfEventModelList.size (); ii++) {
    os << "\n  " << *(d.fEndOfEventModelList[ii]);
  }

  os << "\n  Extent or bounding box: " << d.fExtent;

  os << "\n  Standard target point:  " << d.fStandardTargetPoint;

  return os;
}

G4bool G4Scene::operator != (const G4Scene& s) const {
  if (
      (fRunDurationModelList.size () !=
       s.fRunDurationModelList.size ())              ||
      (fExtent               != s.fExtent)              ||
      !(fStandardTargetPoint == s.fStandardTargetPoint)
      ) return true;

  for (int i = 0; i < fRunDurationModelList.size (); i++) {
    if (fRunDurationModelList[i] != s.fRunDurationModelList[i])
      return true;
  }

  return false;
}
