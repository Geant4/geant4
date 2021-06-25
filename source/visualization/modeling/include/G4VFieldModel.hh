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
// Michael Kelsey  31 January 2019
//
// Class Description:
//
// Abstract base class to implement drawing vector field geometries
// (e.g., electric, magnetic or gravity).  Implementation extracted
// from G4MagneticFieldModel, with field-value access left pure
// virtual for implementation by base classes.

#ifndef G4VFIELDMODEL_HH
#define G4VFIELDMODEL_HH

#include "G4VModel.hh"
#include "G4Point3D.hh"
#include "G4PhysicalVolumesSearchScene.hh"

#include <vector>

class G4Field;

class G4VFieldModel: public G4VModel {

public: // With description
  
  enum Representation {fullArrow, lightArrow};

  G4VFieldModel
  (const G4String& typeOfField, const G4String& symbol="",
   const G4VisExtent& extentForField = G4VisExtent(),
   const std::vector<G4PhysicalVolumesSearchScene::Findings>& pvFindings
   = std::vector<G4PhysicalVolumesSearchScene::Findings>(),
   G4int nDataPointsPerHalfScene = 10,
   Representation representation = Representation::fullArrow,
   G4int arrow3DLineSegmentsPerCircle = 6);
  // typeOfField is "Electric" or "Magnetic" etc.
  // symbol is "E" or "B" etc.

  virtual ~G4VFieldModel();

  virtual void DescribeYourselfTo(G4VGraphicsScene& sceneHandler);
  // The main task of a model is to describe itself to the graphics scene.
  // Note: It is in this function that the extent for drawing the filed must
  // be calcualted. If fExtentForField is null, pick up the extent from
  // the sceneHandler.

protected:

  // Subclasses MUST implement this for their particular kind of field
  virtual void GetFieldAtLocation(const G4Field* field,
				  const G4Point3D& position, G4double time,
				  G4Point3D& result) const = 0;
  // The appropriate output from GetFieldValue should be filled into result.
  // If (field==0), the function should do nothing; returning without error.

private:

  // Private copy contructor and assignment to forbid use...
  G4VFieldModel(const G4VFieldModel&);
  G4VFieldModel& operator=(const G4VFieldModel&);

  G4VisExtent fExtentForField;
  // If null, get extent from scene handler in DescribeYourselfTo.

  std::vector<G4PhysicalVolumesSearchScene::Findings> fPVFindings;
  // If empty, use fExtentForField alone for sampling and drawing.
  // If non-empty, use fExtentForField alone for sampling, but only
  // draw if sampling point is in the specified physical volume(s).

  G4int fNDataPointsPerMaxHalfExtent;
  // No. of data points sampled per maximum half extent.
  // Note that total number of sampling points can be as high as
  // (2*n+1)^3, which can get very big. However, fields are usually
  // confined to only parts of the scene, so this may not be a problem.
  // Sampling can be further limited with fExtentForField and/or fPVFindings.

  Representation fRepresentation;       // Big arrows or just lines
  G4int fArrow3DLineSegmentsPerCircle;
  G4String fTypeOfField;                // "Electric" or "Magnetic" etc.
  G4String fArrowPrefix;                // For attaching text label to arrows
};

#endif
