// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorScene.hh,v 1.1 1999-01-07 16:15:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// J Kallenbach  27th Aug 1996
// OpenInventor scene - creates OpenInventor Display lists.
// 20 dec 1996  jck  Add HEPVis primitives - trd, box, etc.

#ifndef G4OPENINVENTORSCENE_HH
#define G4OPENINVENTORSCENE_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include <rw/tvhdict.h>
#include <rw/tphdict.h>

#include "G4VScene.hh"

class G4OpenInventor;
class SoSeparator;

// Base class for various OpenInventorScene classes.
class G4OpenInventorScene: public G4VScene {

friend class G4OpenInventorView;

public:
  typedef unsigned long      G4VSolidPointer;
  typedef const G4VPhysicalVolume* G4VPhysicalVolumePointer;
  typedef SoSeparator*       SoSeparatorPointer;
  G4OpenInventorScene (G4OpenInventor& system, const G4String& name = "");
  ~G4OpenInventorScene ();
  void AddPrimitive (const G4Polyline& line);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Polyhedron& p);
  void AddPrimitive (const G4NURBS& nurb);
  void AddPrimitive (const G4Polymarker&);
  
  //
  // Primitives for use of HEPVis
  //
  void AddThis (const G4Box  & Box);
  void AddThis (const G4Tubs & Tubs);
  void AddThis (const G4Cons & Cons);
  void AddThis (const G4Trd  & Trd);
  void AddThis (const G4Trap & Trap);

  // Base class callbacks defined in G4OpenInventorScene.icc
  void AddThis (const G4Sphere&);
  void AddThis (const G4Para&);
  void AddThis (const G4Torus&);
  void AddThis (const G4VSolid&);

  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void EndModeling ();
  static G4int GetSceneCount ();
  void PreAddThis (const G4Transform3D& objectTransformation,
		   const G4VisAttributes& visAttribs);

private:
  void 		ClearStore ();
  void 		RequestPrimitives (const G4VSolid& solid);
  static 	G4int    fSceneIdCount;   // static counter for OpenInventor scenes.
  static 	G4int    fSceneCount;     // No. of extanct scenes.
  G4double  	GetMarkerSize    ( const G4VMarker&  mark ) ;

  //
  // Stop-gap solution of structure re-use.
  // A proper implementation would use geometry hierarchy.
  //
  RWTValHashDictionary<G4VSolidPointer, G4int> fSolidDictionary;
  RWTValHashDictionary<G4VPhysicalVolumePointer, SoSeparatorPointer> SeparatorSet;
  SoSeparator *root;
  SoSeparator *currentSeparator;
};

#include "G4OpenInventorScene.icc"

#endif

#endif
