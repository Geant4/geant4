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
//
// 
// Michael Kelsey  31st January 2019 -- Adapted from new G4MagneticFieldModel
//
// Class Description:
//
// Model that knows how to draw the electric field.

#include "G4ElectricFieldModel.hh"
#include "G4Field.hh"
#include "G4Point3D.hh"


// Return electric field vector for display

void G4ElectricFieldModel::
GetFieldAtLocation(const G4Field* field, const G4Point3D& position,
		   G4double time, G4Point3D& result) const {
  if (!field) return;			// No action if no field

  G4double xyzt[4] = { position.x(), position.y(), position.z(), time };
  G4double BEvals[6] = {0.};		// Field returns {Bx,By,Bz,Ex,Ey,Ez}
  field->GetFieldValue(xyzt, BEvals);

  result.set(BEvals[3], BEvals[4], BEvals[5]);
  return;
}
