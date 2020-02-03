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

#ifndef G4ELECTRICFIELDMODEL_HH
#define G4ELECTRICFIELDMODEL_HH

#include "G4VFieldModel.hh"


class G4ElectricFieldModel: public G4VFieldModel {
public: 		// With description

  // Constructor just passes through to base
  G4ElectricFieldModel
  (G4int nDataPointsPerHalfExtent = 3,
   Representation representation = Representation::fullArrow,
   G4int arrow3DLineSegmentsPerCircle = 6,
   const G4VisExtent& extentForField = G4VisExtent(),
   const std::vector<G4PhysicalVolumesSearchScene::Findings>& pvFindings
   = std::vector<G4PhysicalVolumesSearchScene::Findings>())
  : G4VFieldModel
  ("Electric","E", extentForField, pvFindings,
   nDataPointsPerHalfExtent, representation, arrow3DLineSegmentsPerCircle)
  {}

  virtual ~G4ElectricFieldModel() {;}

protected:
  virtual void GetFieldAtLocation(const G4Field* field,
				  const G4Point3D& position, G4double time,
				  G4Point3D& result) const;
  // The appropriate output from GetFieldValue should be filled into result.
  // If (field==0), the function should do nothing; returning without error.

private:
  // Private copy contructor and assignment to forbid use...
  G4ElectricFieldModel(const G4ElectricFieldModel&);
  G4ElectricFieldModel& operator=(const G4ElectricFieldModel&);
};

#endif
