// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VModel.hh,v 1.3 1999-11-11 15:38:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  31st December 1997.
// Base class for models.

#ifndef G4VMODEL_HH
#define G4VMODEL_HH

#include "globals.hh"
#include "G4VisExtent.hh"
#include "G4Transform3D.hh"

class G4VGraphicsScene;
class G4ModelingParameters;

class G4VModel {

  friend ostream& operator << (ostream& os, const G4VModel&);

public:

  G4VModel
  (const G4Transform3D& modelTransformation = G4Transform3D::Identity,
   const G4ModelingParameters* = 0);
   
  virtual ~G4VModel ();

  // For G4RWTPtrOrderedVector...
  G4bool operator == (const G4VModel&) const;

  virtual void DescribeYourselfTo (G4VGraphicsScene&) = 0;
  // The main task of a model is to describe itself to the scene.

  virtual G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

  virtual G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  virtual G4bool Validate ();
  // Validate, but allow internal changes (hence non-const function).

  //////////////////////////////////////////////////////////
  // Access functions...

  const G4ModelingParameters* GetModelingParameters () const;

  const G4String& GetGlobalTag () const;
  // A tag which does not change and lasts the life of the model.
  // Define protected data member in derived class constructor.

  const G4String& GetGlobalDescription () const;
  // A description which does not change and lasts the life of the model.
  // Define protected data member in derived class constructor.

  const G4VisExtent& GetExtent () const;
  // Extent of visible objects in local coordinate system.
  // Define protected data member in derived class constructor.

  const G4Transform3D& GetTransformation () const;
  // Model transformation, i.e., position and orientation of model in world.

  void SetModelingParameters (const G4ModelingParameters*);

protected:

  G4String                    fGlobalTag;
  G4String                    fGlobalDescription;
  const G4ModelingParameters* fpMP;
  G4VisExtent                 fExtent;
  G4Transform3D               fTransform;           

};

#include "G4VModel.icc"

#endif
