// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoriesModel.hh,v 1.1 1999-01-07 16:15:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 trajectories.

#ifndef G4TRAJECTORIESMODEL_HH
#define G4TRAJECTORIESMODEL_HH

#include "G4VModel.hh"

class G4TrajectoriesModel: public G4VModel {

public:

  G4TrajectoriesModel ();
   
  virtual ~G4TrajectoriesModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the scene.

  virtual G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

  virtual G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  virtual G4bool Validate ();
  // Validate, but allow internal changes (hence non-const function).

};

#include "G4TrajectoriesModel.icc"

#endif
