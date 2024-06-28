//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// 
// J Kallenbach  27th Aug 1996
// OpenInventor scene handler - creates OpenInventor Display lists.
// 20 dec 1996  jck  Add HEPVis primitives - trd, box, etc.

#ifndef G4OPENINVENTORSCENEHANDLER_HH
#define G4OPENINVENTORSCENEHANDLER_HH

// Inheritance :
#include "G4VSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"

#include <map>

class G4OpenInventor;
class SoSeparator;
class Geant4_SoStyleCache;
class G4VisAttributes;

// Base class for various OpenInventorScene classes.
class G4OpenInventorSceneHandler: public G4VSceneHandler {

friend class G4OpenInventorViewer;

public:

  G4OpenInventorSceneHandler (G4OpenInventor& system, const G4String& name = "");
  virtual ~G4OpenInventorSceneHandler ();

  using G4VSceneHandler::AddPrimitive;
  void AddPrimitive (const G4Polyline& line);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Polyhedron& p);
  void AddPrimitive (const G4Polymarker&);

  using G4VSceneHandler::AddCompound;
  void AddCompound (const G4Mesh&);

  ///////////////////////////////////////////////////////////////
  // Other inherited functions.
  void 		ClearStore ();
  void 		ClearTransientStore ();
  
  //
  // Primitives for use of HEPVis
  //
  void PreAddSolid (const G4Transform3D& objectTransformation,
		    const G4VisAttributes& visAttribs);
  void BeginPrimitives (const G4Transform3D& objectTransformation);

private:

  static G4int fSceneIdCount;   // static counter for OpenInventor scenes.
  enum G4OIMarker {G4OICircle, G4OISquare};
  void AddCircleSquare (G4OIMarker markerType, const G4VMarker&);
  void GeneratePrerequisites();
  void AddProperties(const G4VisAttributes*);
  // AddTransform takes fObjectTransformation and "adds" a translation.
  void AddTransform(const G4Point3D& translation = G4Point3D());
  std::map <G4LogicalVolume*, SoSeparator*,
    std::less <G4LogicalVolume*> > fSeparatorMap;
  SoSeparator* fRoot;
  SoSeparator* fDetectorRoot;
  SoSeparator* fTransientRoot;
  SoSeparator* fCurrentSeparator;
  G4bool fModelingSolid;
  G4bool fReducedWireFrame;
  Geant4_SoStyleCache* fStyleCache;
  bool fPreviewAndFull;
};

#endif
