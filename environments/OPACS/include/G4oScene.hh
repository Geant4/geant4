// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4oScene.hh,v 1.2 1999/04/16 10:03:27 barrand Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
/* +---------------------- Copyright notice -------------------------------+ */
/* | Copyright (C) 1995, Guy Barrand, LAL Orsay, (barrand@lal.in2p3.fr)    | */
/* |   Permission to use, copy, modify, and distribute this software       | */
/* |   and its documentation for any purpose and without fee is hereby     | */
/* |   granted, provided that the above copyright notice appear in all     | */
/* |   copies and that both that copyright notice and this permission      | */
/* |   notice appear in supporting documentation.  This software is        | */
/* |   provided "as is" without express or implied warranty.               | */
/* +---------------------- Copyright notice -------------------------------+ */

#ifndef G4OSCENE_HH
#define G4OSCENE_HH

#include <G4VGraphicsScene.hh>

class G4PhysicalVolumeModel;
class G4VPhysicalVolume;
class G4LogicalVolume;

#include <ONode.h>

class G4oScene : public G4VGraphicsScene {
public:
  //Inherited
  void AddThis (const G4Box&);
  void AddThis (const G4Cons&);
  void AddThis (const G4Tubs&);
  void AddThis (const G4Trd&);
  void AddThis (const G4Trap&);
  void AddThis (const G4Sphere&);
  void AddThis (const G4Para&);
  void AddThis (const G4Torus&);
  void AddThis (const G4Polycone&);
  void AddThis (const G4Polyhedra&);
  void AddThis (const G4VSolid& solid);  // For solids not above.

  void PreAddThis        (const G4Transform3D& objectTransformation,
			  const G4VisAttributes& visAttribs);
  void PostAddThis       ();

  void EstablishSpecials (G4PhysicalVolumeModel&);

  //Local
  G4oScene ();
  ~G4oScene ();
  void      SetNodeName (char*);
  ONode     GetNode    ();

private:
  static ONode         node;
  const G4Transform3D* fpAccTransf;    // Accumulated transformation.
  static char*         nodeName;
  void AddPolyhedron   (G4Polyhedron*);
  G4int                fCurrentDepth; // Current depth of geom. hierarchy.
  G4VPhysicalVolume*   fpCurrentPV;   // Current physical volume.
  G4LogicalVolume*     fpCurrentLV;   // Current logical volume.
};
#endif
