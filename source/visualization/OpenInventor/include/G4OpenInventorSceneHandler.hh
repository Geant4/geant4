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
// $Id: G4OpenInventorSceneHandler.hh,v 1.23 2004/11/24 14:59:38 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 
// J Kallenbach  27th Aug 1996
// OpenInventor scene handler - creates OpenInventor Display lists.
// 20 dec 1996  jck  Add HEPVis primitives - trd, box, etc.

#ifndef G4OPENINVENTORSCENEHANDLER_HH
#define G4OPENINVENTORSCENEHANDLER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

// Inheritance :
#include "G4VSceneHandler.hh"

#include <map>

class G4OpenInventor;
class SoSeparator;
class Geant4_SoStyleCache;

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
  ////////////////////////////////////////////////////////////////
  // Explicitly invoke base class methods to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive (const G4Scale& scale) {
    G4VSceneHandler::AddPrimitive (scale);
  }
  
  //
  // Primitives for use of HEPVis
  //
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
  G4double  	GetMarkerSize    ( const G4VMarker&  mark ) ;
private:
  static G4int fSceneIdCount;   // static counter for OpenInventor scenes.
  static G4int fSceneCount;
private:
  //
  // Stop-gap solution of structure re-use.
  // A proper implementation would use geometry hierarchy.
  //
  std::map <const G4LogicalVolume*, SoSeparator*,
    std::less <const G4LogicalVolume*> > fSeparatorMap;
  SoSeparator* fRoot;
  SoSeparator* fDetectorRoot;
  SoSeparator* fTransientRoot;
  SoSeparator* fCurrentSeparator;
  G4bool fModelingSolid;
  G4bool fReducedWireFrame;
  Geant4_SoStyleCache* fStyleCache;
  bool fPreviewAndFull;
};

#include "G4OpenInventorSceneHandler.icc"

#endif

#endif
