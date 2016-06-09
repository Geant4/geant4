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
// $Id: G4VModel.cc,v 1.11 2005/11/30 16:06:18 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// John Allison  31st December 1997.
// Base class for models.

#include "G4VModel.hh"

#include "G4RotationMatrix.hh"
#include "G4ModelingParameters.hh"

G4VModel::G4VModel (const G4Transform3D& modelTransformation,
		    const G4ModelingParameters* pMP):
  fGlobalTag ("Empty"),
  fGlobalDescription ("Empty"),
  fTransform (modelTransformation),
  fpMP (pMP)
{}

G4VModel::~G4VModel () {}

G4String G4VModel::GetCurrentTag () const {
  return G4String
    ("WARNING: GetCurrentTag() not implemented by concrete class."
     "\n  Global tag: " + fGlobalTag);
}

G4String G4VModel::GetCurrentDescription () const {
  return G4String
    ("WARNING: GetCurrentDescription() not implemented by concrete class."
     "\n  Global description: " + fGlobalDescription);
}

G4bool G4VModel::Validate (G4bool) {
  return true;
}

std::ostream& operator << (std::ostream& os, const G4VModel& m) {
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
  G4RotationMatrix rotation = m.fTransform.getRotation ();
  os << rotation.thetaX() << ", "
     << rotation.phiX() << ", "
     << rotation.thetaY() << ", "
     << rotation.phiY() << ", "
     << rotation.thetaZ() << ", "
     << rotation.phiZ();
  os << "\n    Translation: " << m.fTransform.getTranslation ();
  return os;
}
