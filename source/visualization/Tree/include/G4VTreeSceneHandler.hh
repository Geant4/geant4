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
// $Id: G4VTreeSceneHandler.hh,v 1.6 2001-08-24 20:41:29 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// In the derived class, override G4VScenehandler::RequestPrimitives
// to implement dump of the "leaves" of the geometry heirachy.

#ifndef G4VTREESCENEHANDLER_HH
#define G4VTREESCENEHANDLER_HH

#include "G4VSceneHandler.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4ModelingParameters;

class G4VTreeSceneHandler: public G4VSceneHandler {

public:
  G4VTreeSceneHandler(G4VGraphicsSystem& system,
		      const G4String& name);
  virtual ~G4VTreeSceneHandler ();
  void PreAddThis (const G4Transform3D& objectTransformation,
		   const G4VisAttributes&);
  void PostAddThis ();
  void EstablishSpecials (G4PhysicalVolumeModel&);
  G4int                GetFoundDepth          () const;
  G4VPhysicalVolume*   GetFoundVolume         () const;
  const G4Transform3D& GetFoundTransformation () const;

  ////////////////////////////////////////////////////////////////
  // Functions not used but required by the abstract interface.

  virtual void BeginPrimitives (const G4Transform3D& objectTransformation) {}
  virtual void EndPrimitives () {}
  virtual void AddPrimitive (const G4Polyline&)   {}
  virtual void AddPrimitive (const G4Text&)       {}
  virtual void AddPrimitive (const G4Circle&)     {}
  virtual void AddPrimitive (const G4Square&)     {}
  virtual void AddPrimitive (const G4Polyhedron&) {}
  virtual void AddPrimitive (const G4NURBS&)      {}

  virtual void BeginModeling();
  virtual void EndModeling();

  static G4int GetSceneCount();

protected:
  // In the derived class, override G4VScenehandler::RequestPrimitives
  // to implement dump of the "leaves" of the geometry heirachy.
  static G4int         fSceneIdCount;  // Counter for Tree scene handlers.
  static G4int         fSceneCount;    // No. of extanct scene handlers.
  G4int                fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume*   fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*     fpCurrentLV;    // Current logical volume.
  const G4Transform3D* fpCurrentObjectTransformation;
  const G4ModelingParameters* fpOriginalMP;  // Keeps pointer to original.
  G4ModelingParameters* fpNonCullingMP;      // For temporary non-culling.
};

#include "G4VTreeSceneHandler.icc"

#endif
