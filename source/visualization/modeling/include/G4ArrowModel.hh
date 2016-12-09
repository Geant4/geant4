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
// $Id: G4ArrowModel.hh 100807 2016-11-02 15:00:41Z gcosmo $
//
// 
// John Allison  15th July 2012
//
// Class Description:
//
// Model that knows how to draw an arrow.
//
// For access to base class information, e.g., modeling parameters,
// use GetModelingParameters() inherited from G4VModel.  See Class
// Description of the base class G4VModel.

#ifndef G4ARROWMODEL_HH
#define G4ARROWMODEL_HH

#include "G4VModel.hh"

class G4Colour;
class G4Polyhedron;

class G4ArrowModel: public G4VModel {

public: // With description

  G4ArrowModel(G4double x1, G4double y1, G4double z1,
	       G4double x2, G4double y2, G4double z2,
	       G4double width, const G4Colour& colour,
	       const G4String& description = "",
               G4int lineSegmentsPerCircle = 6);
  virtual ~G4ArrowModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene.

private:

  // Private copy contructor and assignment to forbid use...
  G4ArrowModel (const G4ArrowModel&);
  G4ArrowModel& operator = (const G4ArrowModel&);

  G4Polyhedron* fpShaftPolyhedron;
  G4Polyhedron* fpHeadPolyhedron;
};

#endif
