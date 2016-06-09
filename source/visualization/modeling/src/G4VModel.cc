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
// $Id: G4VModel.cc,v 1.14 2006-11-01 10:33:00 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  // Override in concrete class if concept of "current" is meaningful.
  return fGlobalTag;
}

G4String G4VModel::GetCurrentDescription () const {
  // Override in concrete class if concept of "current" is meaningful.
  return fGlobalDescription;
}

G4bool G4VModel::Validate (G4bool) {
  return true;
}

std::ostream& operator << (std::ostream& os, const G4VModel& m) {
  os << m.fGlobalDescription;
  os << "\n  Modeling parameters:";
  const G4ModelingParameters* mp = m.fpMP;
  if (mp) os << "\n  " << *mp;
  else os << " none.";
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
