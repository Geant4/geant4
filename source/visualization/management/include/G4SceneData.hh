// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneData.hh,v 1.1 1999-01-07 16:15:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Scene data  John Allison  19th July 1996.
//
// Defines the scene in terms of GEANT4 objects.  All other parameters are
// in G4ViewParamaters.
//
// Note that although the C++ name is G4SceneData, it is often simply
// referred to as the scene.

#ifndef G4SCENEDATA_HH
#define G4SCENEDATA_HH

#include "globals.hh"
#include "G4ios.hh"

class G4VPhysicalVolume;

#include "G4VisExtent.hh"
#include "G4Point3D.hh"
#include "G4VModel.hh"
#include <rw/tpordvec.h>

class G4SceneData {

  friend ostream& operator << (ostream& os, const G4SceneData& d);
  friend G4bool operator != (const G4SceneData& d1,
			     const G4SceneData& d2);
public:

  enum {UNLIMITED = -1};

  G4SceneData (const G4String& name = "invalid scene");
  ~G4SceneData ();

  // Makes use of default (compiler generated) copy constructor and
  // assignment operator.

  // For RWTPtrOrderedVector...
  G4bool operator == (const G4SceneData& sd) const;

  //////////////////////////////////////////////////////
  // Get functions...

  const G4String& GetName () const;

  G4bool IsEmpty () const;

  const RWTPtrOrderedVector <G4VModel>& GetRunDurationModelList () const;
  // Contains models which are expected to last for the duration of
  // the run, for example geometry volumes.

  const RWTPtrOrderedVector <G4VModel>& GetEndOfEventModelList () const;
  // Contains models which are described at the end of event when the
  // scene is current.

  const G4VisExtent& GetExtent () const;
  // Overall extent of all objects in the scene.

  const G4Point3D& GetStandardTargetPoint () const;
  // Usually centre of extent.  See G4ViewParameters for definition.

  //////////////////////////////////////////////
  // Add, Set, Clear functions...

  void AddRunDurationModel (G4VModel*);
  // Adds models of type which are expected to last for the duration
  // of the run, for example geometry volumes.

  void AddEndOfEventModel (G4VModel*);
  // Adds models of type which are described at the end of event when
  // the scene is current.

  RWTPtrOrderedVector <G4VModel>& SetRunDurationModelList ();
  // Allows you to change the model list - do with care!

  RWTPtrOrderedVector <G4VModel>& SetEndOfEventModelList ();
  // Allows you to change the model list - do with care!

  void Clear ();
  // Clears and destroys models in all lists.

private:
  G4String fName;
  RWTPtrOrderedVector <G4VModel> fRunDurationModelList;
  RWTPtrOrderedVector <G4VModel> fEndOfEventModelList;
  G4VisExtent fExtent;
  G4Point3D   fStandardTargetPoint;
};

#include "G4SceneData.icc"

#endif
