// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AxesModel.hh,v 1.1 2001-04-11 13:39:30 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  3rd April 2001
//
// Class Description:
//
// Model which knows how to draw axes.
//
// For access to base class information, e.g., modeling parameters,
// use GetModelingParameters() inherited from G4VModel.  See Class
// Description of the base class G4VModel.

#ifndef G4AXESMODEL_HH
#define G4AXESMODEL_HH

#include "G4VModel.hh"

class G4AxesModel: public G4VModel {

public: // With description

  G4AxesModel (G4double, G4double, G4double, G4double);
   
  virtual ~G4AxesModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the scene.

  virtual G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  virtual G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

  virtual G4bool Validate ();
  // Validate, but allow internal changes (hence non-const function).

private:

  // Private copy contructor and assignmen to forbid uset...
  G4AxesModel (const G4AxesModel&);
  G4AxesModel& operator = (const G4AxesModel&);

  G4double fX0, fY0, fZ0, fLength;

};

#include "G4AxesModel.icc"

#endif
