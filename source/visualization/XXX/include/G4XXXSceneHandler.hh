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
// $Id: G4XXXSceneHandler.hh,v 1.16 2005/06/07 16:46:33 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.

#ifndef G4XXXSCENEHANDLER_HH
#define G4XXXSCENEHANDLER_HH

#define G4XXXDEBUG  // Comment this out to suppress debug code.

#include "G4VSceneHandler.hh"

class G4XXXSceneHandler: public G4VSceneHandler {

public:
  G4XXXSceneHandler(G4VGraphicsSystem& system,
		      const G4String& name);
  virtual ~G4XXXSceneHandler();

  ////////////////////////////////////////////////////////////////
  // No need to implement these, but if you do...
  void AddSolid(const G4Box&);
  void AddSolid(const G4Cons&);
  void AddSolid(const G4Tubs&);
  void AddSolid(const G4Trd&);
  void AddSolid(const G4Trap&);
  void AddSolid(const G4Sphere&);
  void AddSolid(const G4Para&);
  void AddSolid(const G4Torus&);
  void AddSolid(const G4Polycone&);
  void AddSolid(const G4Polyhedra&);
  void AddSolid(const G4VSolid&);
  void AddCompound(const G4VTrajectory&);
  void AddCompound(const G4VHit&);
  void PreAddSolid(const G4Transform3D& objectTransformation,
		   const G4VisAttributes&);
  // void PostAddSolid();

  ////////////////////////////////////////////////////////////////
  // Required implementation of pure virtual functions...

  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);
  void AddPrimitive(const G4NURBS&);

  ////////////////////////////////////////////////////////////////
  // Further optional AddPrimtive methods.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive(const G4Polymarker& polymarker) {
    G4VSceneHandler::AddPrimitive (polymarker);
  }
  void AddPrimitive(const G4Scale& scale) {
    G4VSceneHandler::AddPrimitive (scale);
  }

  ////////////////////////////////////////////////////////////////
  // Further optional AddPrimtive methods...

  // void AddPrimitive(const G4Scale&);
  // void AddPrimitive(const G4Polymarker&);

  ////////////////////////////////////////////////////////////////
  // Further optional virtual functions...

  // void BeginPrimitives(const G4Transform3D& objectTransformation);
  // void EndPrimitives();

  // void BeginModeling();
  // void EndModeling();

  //////////////////////////////////////////////////////////////
  // Administration functions.

  //void ClearStore ();
  void ClearTransientStore ();

protected:
  static G4int         fSceneIdCount;  // Counter for XXX scene handlers.

private:
#ifdef G4XXXDEBUG
  void PrintThings();
#endif
};

#endif
