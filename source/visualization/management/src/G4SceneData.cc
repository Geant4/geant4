// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneData.cc,v 1.1 1999-01-07 16:15:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Scene data  John Allison  19th July 1996.

#include "G4SceneData.hh"

#include "G4Vector3D.hh"
#include "G4BoundingSphereScene.hh"

G4SceneData::G4SceneData (const G4String& name):
  fName (name)
{} // Note all data members have default initial values.

G4SceneData::~G4SceneData () {}

void G4SceneData::AddRunDurationModel (G4VModel* pModel) {
  fRunDurationModelList.append (pModel);
  G4int nModels = fRunDurationModelList.entries ();
  G4BoundingSphereScene boundingSphereScene;
  for (int i = 0; i < nModels; i++) {
    const G4VisExtent& thisExtent =
      fRunDurationModelList[i] -> GetExtent ();
    G4Point3D thisCentre = thisExtent.GetExtentCentre ();
    G4double thisRadius = thisExtent.GetExtentRadius ();

    //thisCentre.transform (fRunDurationModelList[i] -> GetTransformation ());
    // To by-pass temporary CLHEP non-const problem...
    G4Transform3D modelTransformation =
      fRunDurationModelList[i] -> GetTransformation ();
    thisCentre.transform (modelTransformation);

    boundingSphereScene.AccrueBoundingSphere (thisCentre, thisRadius);
  }
  fExtent = boundingSphereScene.GetBoundingSphereExtent ();
  fStandardTargetPoint = fExtent.GetExtentCentre ();
}
  
void G4SceneData::Clear () {
  fRunDurationModelList.clearAndDestroy ();
  fEndOfEventModelList.clearAndDestroy ();
}

ostream& operator << (ostream& os, const G4SceneData& d) {

  os << "Scene data:";

  os << "\n  Run-duration model list:";
  for (int i = 0; i < d.fRunDurationModelList.entries (); i++) {
    os << "\n  " << *(d.fRunDurationModelList[i]);
  }

  os << "\n  End-of-event model list:";
  for (int ii = 0; ii < d.fEndOfEventModelList.entries (); ii++) {
    os << "\n  " << *(d.fEndOfEventModelList[ii]);
  }

  os << "\n  Extent or bounding box: " << d.fExtent;

  os << "\n  Standard target point:  " << d.fStandardTargetPoint;

  return os;
}

G4bool operator != (const G4SceneData& d1,
		    const G4SceneData& d2) {
  if (
      (d1.fRunDurationModelList.entries () !=
       d2.fRunDurationModelList.entries ())                 ||
      (d1.fExtent               != d2.fExtent)              ||
      !(d1.fStandardTargetPoint == d2.fStandardTargetPoint)
      ) return true;

  for (int i = 0; i < d1.fRunDurationModelList.entries (); i++) {
    if (d1.fRunDurationModelList[i] != d2.fRunDurationModelList[i])
      return true;
  }

  return false;
}
