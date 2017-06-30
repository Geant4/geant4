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
// $Id: G4MagneticFieldModel.hh 101958 2016-12-12 08:04:35Z gcosmo $
//
// 
// John Allison  17th August 2013
//
// Class Description:
//
// Model that knows how to draw the magnetic field.

#ifndef G4MAGNETICFIELDMODEL_HH
#define G4MAGNETICFIELDMODEL_HH

#include "G4VModel.hh"

class G4Colour;
class G4Polyhedron;

class G4MagneticFieldModel: public G4VModel {

public: // With description

  enum Representation {fullArrow, lightArrow};

  G4MagneticFieldModel
  (G4int nDataPointsPerHalfScene = 10,
   Representation representation = Representation::fullArrow,
   G4int arrow3DLineSegmentsPerCircle = 6);
  virtual ~G4MagneticFieldModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene.

private:

  // Private copy contructor and assignment to forbid use...
  G4MagneticFieldModel (const G4MagneticFieldModel&);
  G4MagneticFieldModel& operator = (const G4MagneticFieldModel&);

  // No. of data points sampled per maximum half scene extent.
  // Note that total number of data poinrs sampled can be as high as
  // (2*n+1)^3, which can get very big very soon.
  G4int fNDataPointsPerMaxHalfScene;

  Representation fRepresentation;

  G4int fArrow3DLineSegmentsPerCircle;
};

#endif
