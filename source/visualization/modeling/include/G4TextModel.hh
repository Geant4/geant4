// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TextModel.hh,v 1.1 2001-04-11 13:39:32 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  3rd April 2001
//
// Class Description:
//
// Model which knows how to draw text.
//
// For access to base class information, e.g., modeling parameters,
// use GetModelingParameters() inherited from G4VModel.  See Class
// Description of the base class G4VModel.

#ifndef G4TEXTMODEL_HH
#define G4TEXTMODEL_HH

#include "G4VModel.hh"
#include "G4Text.hh"

class G4TextModel: public G4VModel {

public: // With description

  G4TextModel (const G4Text&);
   
  virtual ~G4TextModel ();

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
  G4TextModel (const G4TextModel&);
  G4TextModel& operator = (const G4TextModel&);

  G4Text fText;

};

#include "G4TextModel.icc"

#endif
