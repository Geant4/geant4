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
// $Id:$
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4UMultiUnion
//
// Class description:
//
//   Wrapper class for G4UMultiUnion to make use of it from USolids module.

// History:
// 30.10.13 G.Cosmo, T.Nikitina, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UMULTIUNION_HH
#define G4UMULTIUNION_HH

#include <CLHEP/Vector/Rotation.h>

#include "G4USolid.hh"
#include "UMultiUnion.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "HepPolyhedronProcessor.h"

class G4UMultiUnion : public G4USolid 
{
  public:  // with description

    G4UMultiUnion(const G4String& name);
   ~G4UMultiUnion();

    inline UMultiUnion* GetShape() const;

    inline void AddNode(G4VSolid& solid, G4Transform3D& trans);
      // Build the multiple union by adding nodes
    inline G4Transform3D* GetTransformation(G4int index) const;
    inline G4VSolid* GetSolid(G4int index) const;
    inline int GetNumberOfSolids()const;
    inline void Voxelize();
  public:  // without description

    G4UMultiUnion(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UMultiUnion( const G4UMultiUnion& source );
    G4UMultiUnion& operator=(const G4UMultiUnion& source);
      // Copy constructor and assignment operator.

    G4Polyhedron* CreatePolyhedron() const;
      // Called by visualization engine
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UMultiUnion* G4UMultiUnion::GetShape() const
{
  return (UMultiUnion*) fShape;
}

inline void G4UMultiUnion::AddNode(G4VSolid& solid, G4Transform3D& trans)
{
  HepGeom::Rotate3D rot;
  HepGeom::Translate3D transl ;
  HepGeom::Scale3D scale;

  trans.getDecomposition(scale,rot,transl); 
  G4ThreeVector pos = transl.getTranslation();
    
  UTransform3D tr;
  tr.fRot[0] = rot.xx(); tr.fRot[1] = rot.xy(); tr.fRot[2] = rot.xz();
  tr.fRot[3] = rot.yx(); tr.fRot[4] = rot.yy(); tr.fRot[5] = rot.yz();
  tr.fRot[6] = rot.zx(); tr.fRot[7] = rot.zy(); tr.fRot[8] = rot.zz();
  tr.fTr = UVector3(pos.x(), pos.y(), pos.z());
 
  GetShape()->AddNode(*(static_cast<G4USolid&>(solid).GetSolid()), tr);
}

inline G4Transform3D* G4UMultiUnion::GetTransformation(G4int index) const
{
  UTransform3D tr = GetShape()->GetTransformation(index);

  G4RotationMatrix
    rot(CLHEP::HepRep3x3(tr.fRot[0], tr.fRot[1], tr.fRot[2],
                         tr.fRot[3], tr.fRot[4], tr.fRot[5],
                         tr.fRot[6], tr.fRot[7], tr.fRot[8]));
  G4ThreeVector transl(tr.fTr.x(), tr.fTr.y(), tr.fTr.z());

  return new G4Transform3D(rot, transl);
}

inline G4VSolid* G4UMultiUnion::GetSolid(G4int index) const
{
  VUSolid* solid = GetShape()->GetSolid(index);
  return new G4USolid(solid->GetName(), solid);
}

inline int  G4UMultiUnion::GetNumberOfSolids()const
{
  return GetShape()->GetNumberOfSolids();
}

inline void G4UMultiUnion::Voxelize()
{
  GetShape()->Voxelize();
}
#endif
