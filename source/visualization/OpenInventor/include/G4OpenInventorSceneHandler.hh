// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorSceneHandler.hh,v 1.9 2000-08-19 18:34:42 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// J Kallenbach  27th Aug 1996
// OpenInventor scene handler - creates OpenInventor Display lists.
// 20 dec 1996  jck  Add HEPVis primitives - trd, box, etc.

#ifndef G4OPENINVENTORSCENEHANDLER_HH
#define G4OPENINVENTORSCENEHANDLER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include "g4std/map"

#include "G4VSceneHandler.hh"

class G4OpenInventor;
class SoSeparator;

// Base class for various OpenInventorScene classes.
class G4OpenInventorSceneHandler: public G4VSceneHandler {

friend class G4OpenInventorViewer;

public:

  G4OpenInventorSceneHandler (G4OpenInventor& system, const G4String& name = "");
  virtual ~G4OpenInventorSceneHandler ();
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

  // Base class callbacks defined in G4OpenInventorSceneHandler.icc
  void AddThis (const G4Sphere&);
  void AddThis (const G4Para&);
  void AddThis (const G4Torus&);
  void AddThis (const G4Polycone&);
  void AddThis (const G4Polyhedra&);
  void AddThis (const G4VSolid&);

  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void EndModeling ();
  static G4int GetSceneCount ();
  void PreAddThis (const G4Transform3D& objectTransformation,
		   const G4VisAttributes& visAttribs);

private:
  void 		ClearStore ();
  void 		ClearTransientStore ();
  void 		RequestPrimitives (const G4VSolid& solid);
  static 	G4int    fSceneIdCount;   // static counter for OpenInventor scenes.
  static 	G4int    fSceneCount;     // No. of extanct scene handlers.
  G4double  	GetMarkerSize    ( const G4VMarker&  mark ) ;

  //
  // Stop-gap solution of structure re-use.
  // A proper implementation would use geometry hierarchy.
  //
  G4std::map <const G4VPhysicalVolume*, SoSeparator*,
    G4std::less <const G4VPhysicalVolume*> > SeparatorMap;
  SoSeparator *root;
  SoSeparator *staticRoot;
  SoSeparator *transientRoot;
  SoSeparator *currentSeparator;
};

#include "G4OpenInventorSceneHandler.icc"

#endif

#endif
