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
// $Id: G4ScaleModel.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// John Allison  21st July 2001.
//
// Class Description:
//
// Model which knows how to draw a scale.
//
// For access to base class information, e.g., modeling parameters,
// use GetModelingParameters() inherited from G4VModel.  See Class
// Description of the base class G4VModel.

#ifndef G4SCALEMODEL_HH
#define G4SCALEMODEL_HH

#include "G4VModel.hh"
#include "G4Scale.hh"

class G4ScaleModel: public G4VModel {

public: // With description

  G4ScaleModel (const G4Scale&);
   
  virtual ~G4ScaleModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene.

private:

  // Private copy contructor and assignment to forbid use...
  G4ScaleModel (const G4ScaleModel&);
  G4ScaleModel& operator = (const G4ScaleModel&);

  G4Scale fScale;
};

#endif
