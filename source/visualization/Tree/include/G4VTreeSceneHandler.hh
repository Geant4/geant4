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
// $Id: G4VTreeSceneHandler.hh 66870 2013-01-14 23:38:59Z adotti $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// In the derived class, override G4VScenehandler::RequestPrimitives
// to implement dump of the geometry hierarchy.

#ifndef G4VTREESCENEHANDLER_HH
#define G4VTREESCENEHANDLER_HH

#include "G4VSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include <vector>
#include <set>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4ModelingParameters;

class G4VTreeSceneHandler: public G4VSceneHandler {

public:
  G4VTreeSceneHandler(G4VGraphicsSystem& system,
		      const G4String& name);
  virtual ~G4VTreeSceneHandler ();
  void PreAddSolid (const G4Transform3D& objectTransformation,
		   const G4VisAttributes&);
  void PostAddSolid ();

  ////////////////////////////////////////////////////////////////
  // Functions not used but required by the abstract interface.

  virtual void AddPrimitive (const G4Polyline&)   {}
  virtual void AddPrimitive (const G4Text&)       {}
  virtual void AddPrimitive (const G4Circle&)     {}
  virtual void AddPrimitive (const G4Square&)     {}
  virtual void AddPrimitive (const G4Polyhedron&) {}
  virtual void AddPrimitive (const G4Polymarker&) {}
  virtual void AddPrimitive (const G4Scale&)      {}

  virtual void BeginModeling();
  virtual void EndModeling();

protected:

  // In the derived class, override G4VScenehandler::RequestPrimitives
  // to implement dump of the geometry hierarchy.
  static G4int         fSceneIdCount;  // Counter for Tree scene handlers.
  const G4Transform3D* fpCurrentObjectTransformation;
  std::set<G4LogicalVolume*> fDrawnLVStore;  // Stores encountered LVs.
};

#include "G4VTreeSceneHandler.icc"

#endif
