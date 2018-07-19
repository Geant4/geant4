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
// $Id: G4XXXSGSceneHandler.hh 99152 2016-09-07 08:04:30Z gcosmo $
//
// 
// John Allison  10th March 2006
// A template for a sophisticated graphics driver with a scene graph.
//?? Lines beginning like this require specialisation for your driver.

#ifndef G4XXXSGSCENEHANDLER_HH
#define G4XXXSGSCENEHANDLER_HH

#include "G4VSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include <iostream>

namespace JA {
  // Ad hoc tree class and utilities.
  
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;

  struct Node {
    Node(PVNodeID pvNodeID = PVNodeID(), G4int index = -1):
      fPVNodeID(pvNodeID), fIndex(index) {}
    PVNodeID fPVNodeID;
    G4int fIndex;
    std::vector<Node*> fDaughters;
  };

  void Insert(const PVNodeID* pvPath, size_t pathLength,
	      G4int index, Node* node);

  void PrintTree(std::ostream&, Node*);

  void Clear(Node*);
}

typedef JA::Node SceneGraph;

class G4XXXSGSceneHandler: public G4VSceneHandler {

  friend class G4XXXSGViewer;

public:
  G4XXXSGSceneHandler(G4VGraphicsSystem& system,
			const G4String& name);
  virtual ~G4XXXSGSceneHandler();

  ////////////////////////////////////////////////////////////////
  // Optional virtual functions...
  void AddSolid(const G4Box&);
  // Further optional AddSolid functions.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddSolid(const G4Cons& cons)
  {G4VSceneHandler::AddSolid(cons);}
  void AddSolid(const G4Tubs& tubs)
  {G4VSceneHandler::AddSolid(tubs);}
  void AddSolid(const G4Trd& trd)
  {G4VSceneHandler::AddSolid(trd);}
  void AddSolid(const G4Trap& trap)
  {G4VSceneHandler::AddSolid(trap);}
  void AddSolid(const G4Sphere& sphere)
  {G4VSceneHandler::AddSolid(sphere);}
  void AddSolid(const G4Para& para)
  {G4VSceneHandler::AddSolid(para);}
  void AddSolid(const G4Torus& torus)
  {G4VSceneHandler::AddSolid(torus);}
  void AddSolid(const G4Polycone& polycone)
  {G4VSceneHandler::AddSolid(polycone);}
  void AddSolid(const G4Polyhedra& polyhedra)
  {G4VSceneHandler::AddSolid(polyhedra);}
  void AddSolid(const G4Orb& orb)
  {G4VSceneHandler::AddSolid(orb);}
  void AddSolid(const G4Ellipsoid& ellipsoid)
  {G4VSceneHandler::AddSolid(ellipsoid);}
  void AddSolid(const G4VSolid& solid)
  {G4VSceneHandler::AddSolid(solid);}
  // More optional functions...
  // void AddCompound(const G4VTrajectory&);
  // void AddCompound(const G4VHit&);
  // void AddCompound(const G4THitsMap<G4double>&);
  // void AddCompound(const G4THitsMap<G4StatDouble>&);
  void PreAddSolid(const G4Transform3D& objectTransformation,
		   const G4VisAttributes&);
  void PostAddSolid();

  ////////////////////////////////////////////////////////////////
  // Required implementation of pure virtual functions...

  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);
  // Further optional AddPrimitive methods.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive(const G4Polymarker& polymarker)
  {G4VSceneHandler::AddPrimitive (polymarker);}
  void AddPrimitive(const G4Scale& scale)
  {G4VSceneHandler::AddPrimitive (scale);}
  // Further related optional virtual functions...
  void BeginPrimitives(const G4Transform3D& objectTransformation);
  void EndPrimitives();

  ////////////////////////////////////////////////////////////////
  // Further optional virtual functions...

  // void BeginModeling();
  // void EndModeling();

  //////////////////////////////////////////////////////////////
  // Administration functions.

  void ClearStore ();
  void ClearTransientStore ();

protected:

  static G4int fSceneIdCount;  // Counter for XXXSG scene handlers.

  // Utility for PreAddSolid and BeginPrimitives.
  void CreateCurrentItem(const G4String&);

  //?? Define the scene graph.  (For emulation, use an ad hoc tree class.)
  SceneGraph fSceneGraph;

private:

#ifdef G4XXXFileDEBUG
  void PrintThings();
#endif

};

#endif
