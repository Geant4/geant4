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
// $Id: G4XXXFileSceneHandler.hh,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  7th March 2006
// A template for a file-writing graphics driver.
//?? Lines beginning like this require specialisation for your driver.

#ifndef G4XXXFileSCENEHANDLER_HH
#define G4XXXFileSCENEHANDLER_HH

#include "G4VSceneHandler.hh"

class G4XXXFileSceneHandler: public G4VSceneHandler {

  friend class G4XXXFileViewer;

public:
  G4XXXFileSceneHandler(G4VGraphicsSystem& system,
			const G4String& name);
  virtual ~G4XXXFileSceneHandler();

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
  void AddSolid(const G4VSolid& solid)
  {G4VSceneHandler::AddSolid(solid);}
  // More optional functions...
  // void AddCompound(const G4VTrajectory&);
  // void AddCompound(const G4VHit&);
  // void PreAddSolid(const G4Transform3D& objectTransformation,
  //	   const G4VisAttributes&);
  // void PostAddSolid();

  ////////////////////////////////////////////////////////////////
  // Required implementation of pure virtual functions...

  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);
  void AddPrimitive(const G4NURBS&);
  // Further optional AddPrimitive methods.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive(const G4Polymarker& polymarker)
  {G4VSceneHandler::AddPrimitive (polymarker);}
  void AddPrimitive(const G4Scale& scale)
  {G4VSceneHandler::AddPrimitive (scale);}
  // Further related optional virtual functions...
  // void BeginPrimitives(const G4Transform3D& objectTransformation);
  // void EndPrimitives();

  ////////////////////////////////////////////////////////////////
  // Further optional virtual functions...

  // void BeginModeling();
  // void EndModeling();

  //////////////////////////////////////////////////////////////
  // Administration functions.

  // void ClearStore ();
  // void ClearTransientStore ();

protected:

  static G4int fSceneIdCount;  // Counter for XXXFile scene handlers.
};

#endif
