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
// $Id: G4ScaleModel.hh,v 1.2 2001-08-24 20:36:18 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  virtual G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  virtual G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

private:

  // Private copy contructor and assignment to forbid use...
  G4ScaleModel (const G4ScaleModel&);
  G4ScaleModel& operator = (const G4ScaleModel&);

  G4Scale fScale;
};

#include "G4ScaleModel.icc"

#endif
