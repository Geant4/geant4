// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VModel.hh,v 1.7 2001-02-23 15:43:33 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  31st December 1997.
//
// Class Description:
//
// G4VModel is a base class for visualization models.  A model is a
// graphics-system-indepedent description of a Geant4 component.
// The key fuctionality of a model is to know how to describe itself
// to a scene handler.  A scene is a collection of models.

#ifndef G4VMODEL_HH
#define G4VMODEL_HH

#include "globals.hh"
#include "G4VisExtent.hh"
#include "G4Transform3D.hh"

class G4VGraphicsScene;
class G4ModelingParameters;

class G4VModel {

public: // With description

  friend G4std::ostream& operator << (G4std::ostream& os, const G4VModel&);

  G4VModel
  (const G4Transform3D& modelTransformation = G4Transform3D::Identity,
   const G4ModelingParameters* = 0);
   
  virtual ~G4VModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&) = 0;
  // The main task of a model is to describe itself to the scene.

  const G4ModelingParameters* GetModelingParameters () const;

  virtual G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  virtual G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

  const G4VisExtent& GetExtent () const;
  // Extent of visible objects in local coordinate system.
  // Define protected data member in derived class constructor.

  const G4String& GetGlobalDescription () const;
  // A description which does not change and lasts the life of the model.
  // Define protected data member in derived class constructor.

  const G4String& GetGlobalTag () const;
  // A tag which does not change and lasts the life of the model.
  // Define protected data member in derived class constructor.

  const G4Transform3D& GetTransformation () const;
  // Model transformation, i.e., position and orientation of model in world.

  void SetModelingParameters (const G4ModelingParameters*);

  virtual G4bool Validate ();
  // Validate, but allow internal changes (hence non-const function).

protected:

  G4String                    fGlobalTag;
  G4String                    fGlobalDescription;
  const G4ModelingParameters* fpMP;
  G4VisExtent                 fExtent;
  G4Transform3D               fTransform;           

private:

  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps CodeWizard happy.
  G4VModel (const G4VModel&);
  G4VModel& operator = (const G4VModel&);
};

#include "G4VModel.icc"

#endif
