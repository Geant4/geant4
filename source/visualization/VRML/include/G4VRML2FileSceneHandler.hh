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
// G4VRML2FileSceneHandler.hh
// Satoshi Tanaka & Yasuhide Sawada

#ifndef G4VRML2FILE_SCENE_HANDLER_HH
#define G4VRML2FILE_SCENE_HANDLER_HH

#include <fstream>

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VSceneHandler.hh"

class G4VRML2File;
class G4VisAttributes;

class G4VRML2FileSceneHandler : public G4VSceneHandler
{
  friend class G4VRML2FileViewer;

  // methods (public)
 public:
  G4VRML2FileSceneHandler(G4VRML2File& system, const G4String& name = "");
  virtual ~G4VRML2FileSceneHandler();

  using G4VSceneHandler::AddSolid;
  void AddSolid(const G4Box&);
  void AddSolid(const G4Cons&);
  void AddSolid(const G4Tubs&);
  void AddSolid(const G4Trd&);
  void AddSolid(const G4Trap&);
  void AddSolid(const G4Sphere&);
  void AddSolid(const G4Para&);
  void AddSolid(const G4Torus&);
  void AddSolid(const G4VSolid&);

  using G4VSceneHandler::AddCompound;

  void BeginPrimitives(const G4Transform3D& objectTransformation);
  void EndPrimitives();

  using G4VSceneHandler::AddPrimitive;
  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Polyhedron&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);

  void ClearTransientStore();  // Used for triggering detector re-drawing.

  void BeginModeling();
  void EndModeling();

  void VRMLBeginModeling();
  void VRMLEndModeling();

  void connectPort();
  void closePort();

  // methods (private)
 private:
  void SendMaterialNode(const G4VisAttributes* pAV);
  void SendMaterialNode();

  void SendLineColor(const G4VisAttributes* pAV);
  void SendMarkerColor(const G4VMarker& mark);
  void SendMarkerWorldPosition(const G4VMarker& mark);

  G4double GetMarkerHalfSize(const G4VMarker& mark);
  void GetMarkerWorldPosition(const G4VMarker& mark, double* pX, double* pY,
                              double* pZ);

  G4bool isConnected() { return fFlagDestOpen; }

  G4bool IsPVPickable() { return fPVPickable; }
  void SetPVPickability(G4bool on_off) { fPVPickable = on_off; }
  G4double SetPVTransparency();
  G4double GetPVTransparency() { return fPVTransparency; }

  // data
 private:
  char fVRMLFileDestDir[256];
  char fVRMLFileName[256];

  G4VRML2File& fSystem;  // Graphics system for this scene.
  G4bool fFlagDestOpen;

  G4int fMaxFileNum;

  G4bool fPVPickable;
  G4double fPVTransparency;

  static G4int fSceneIdCount;

 public:
  std::ofstream fDest;
};

#endif  // G4VRML2FILE_SCENE_HANDLER_HH
