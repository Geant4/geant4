// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VModel.cc,v 1.3.6.1 1999/12/07 20:54:12 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// John Allison  31st December 1997.
// Base class for models.

#include "G4VModel.hh"

#include "G4ModelingParameters.hh"

G4VModel::G4VModel (const G4Transform3D& modelTransformation,
		    const G4ModelingParameters* pMP):
  fTransform (modelTransformation),
  fpMP (pMP)
{
  fGlobalTag = "Default Global Tag";
  fGlobalDescription = "Default Global Description";
}

G4VModel::~G4VModel () {}

ostream& operator << (ostream& os, const G4VModel& m) {
  os << m.fGlobalDescription;
  os << "\n  Modeling parameters:";
  if (m.fpMP) {
    os << "\n  " << *(m.fpMP);
  }
  else {
    os << " none.";
  }    
  os << "\n  Extent: " << m.fExtent;
  os << "\n  Transformation: ";
  os << "\n    Rotation: ";
  HepRotation rotation = m.fTransform.getRotation ();
  os << rotation.thetaX() << ", "
     << rotation.phiX() << ", "
     << rotation.thetaY() << ", "
     << rotation.phiY() << ", "
     << rotation.thetaZ() << ", "
     << rotation.phiZ();
  os << "\n    Translation: " << m.fTransform.getTranslation ();
  return os;
}
