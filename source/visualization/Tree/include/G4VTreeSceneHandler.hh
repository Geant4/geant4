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
// $Id: G4VTreeSceneHandler.hh,v 1.16 2006-03-14 12:37:40 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  virtual void BeginPrimitives (const G4Transform3D&) {}
  virtual void EndPrimitives () {}
  virtual void AddPrimitive (const G4Polyline&)   {}
  virtual void AddPrimitive (const G4Text&)       {}
  virtual void AddPrimitive (const G4Circle&)     {}
  virtual void AddPrimitive (const G4Square&)     {}
  virtual void AddPrimitive (const G4Polyhedron&) {}
  virtual void AddPrimitive (const G4NURBS&)      {}
  virtual void AddPrimitive (const G4Polymarker&) {}
  virtual void AddPrimitive (const G4Scale&)      {}

  virtual void BeginModeling();
  virtual void EndModeling();

  ///////////////////////////////////////////////////////////////
  // Other inherited functions.

  void EstablishSpecials (G4PhysicalVolumeModel&);
  // Used to establish any special relationships between scene and this
  // particular type of model - non-pure, i.e., no requirement to
  // implement.  See G4PhysicalVolumeModel.hh for details.

protected:
  // In the derived class, override G4VScenehandler::RequestPrimitives
  // to implement dump of the geometry hierarchy.
  static G4int         fSceneIdCount;  // Counter for Tree scene handlers.
  const G4Transform3D* fpCurrentObjectTransformation;

  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  PVPath fDrawnPVPath;                 // Path of drawn (non-culled) PVs.
  std::set<PVNodeID> fPVNodeStore;     // Stores encountered PVNodeIDs.
  std::set<G4LogicalVolume*> fDrawnLVStore;  // Stores encountered LVs.
};

#include "G4VTreeSceneHandler.icc"

#endif
