// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HitsModel.hh,v 1.2 1999-01-11 00:48:43 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 hits.

#ifndef G4HITSMODEL_HH
#define G4HITSMODEL_HH

#include "G4VModel.hh"

class G4HitsModel: public G4VModel {

public:

  G4HitsModel ();
   
  virtual ~G4HitsModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the scene.

  virtual G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

  virtual G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  virtual G4bool Validate ();
  // Validate, but allow internal changes (hence non-const function).

};

#include "G4HitsModel.icc"

#endif
