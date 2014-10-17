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
// G4UPolyhedra
//
// Class description:
//
//   Wrapper class for UPolyhedra to make use of it from USolids module.

// History:
// 31.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UPOLYHEDRA_HH
#define G4UPOLYHEDRA_HH

#include "G4USolid.hh"
#include "UPolyhedra.hh"
#include "G4PolyhedraSide.hh"
#include "G4PolyhedraHistorical.hh"
#include "G4Polyhedron.hh"

class G4EnclosingCylinder;
class G4ReduciblePolygon;

class G4UPolyhedra : public G4USolid
{
  public:  // with description

    G4UPolyhedra( const G4String& name, 
                     G4double phiStart,    // initial phi starting angle
                     G4double phiTotal,    // total phi angle
                     G4int numSide,        // number sides
                     G4int numZPlanes,     // number of z planes
               const G4double zPlane[],    // position of z planes
               const G4double rInner[],    // tangent distance to inner surface
               const G4double rOuter[]  ); // tangent distance to outer surface

    G4UPolyhedra( const G4String& name, 
                     G4double phiStart,    // initial phi starting angle
                     G4double phiTotal,    // total phi angle
                     G4int    numSide,     // number sides
                     G4int    numRZ,       // number corners in r,z space
               const G4double r[],         // r coordinate of these corners
               const G4double z[]       ); // z coordinate of these corners

   ~G4UPolyhedra();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline UPolyhedra* GetShape() const;

    inline G4int GetNumSide()     const;
    inline G4double GetStartPhi() const;
    inline G4double GetEndPhi()   const;
    inline G4bool IsOpen()        const;
    inline G4bool IsGeneric()     const;
    inline G4int GetNumRZCorner() const;
    inline G4PolyhedraSideRZ GetCorner( const G4int index ) const;
    inline G4PolyhedraHistorical* GetOriginalParameters() const;
    inline void SetOriginalParameters(G4PolyhedraHistorical* pars);

    inline G4bool Reset();

  public:  // without description

    G4UPolyhedra(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UPolyhedra( const G4UPolyhedra &source );
    G4UPolyhedra &operator=( const G4UPolyhedra &source );
      // Copy constructor and assignment operator.
    G4Polyhedron* CreatePolyhedron() const;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UPolyhedra* G4UPolyhedra::GetShape() const
{
  return (UPolyhedra*) fShape;
}

inline G4int G4UPolyhedra::GetNumSide() const
{
  return GetShape()->GetNumSide();
}
inline G4double G4UPolyhedra::GetStartPhi() const
{
  return GetShape()->GetStartPhi();
}
inline G4double G4UPolyhedra::GetEndPhi() const
{
  return GetShape()->GetEndPhi();
}
inline G4bool G4UPolyhedra::IsOpen() const
{
  return GetShape()->IsOpen();
}
inline G4bool G4UPolyhedra::IsGeneric() const
{
  return GetShape()->IsGeneric();
}
inline G4int G4UPolyhedra::GetNumRZCorner() const
{
  return GetShape()->GetNumRZCorner();
}
inline G4PolyhedraSideRZ G4UPolyhedra::GetCorner(G4int index) const
{
  UPolyhedraSideRZ pside = GetShape()->GetCorner(index);
  G4PolyhedraSideRZ psiderz = { pside.r, pside.z };

  return psiderz;
}
inline G4PolyhedraHistorical* G4UPolyhedra::GetOriginalParameters() const
{
  UPolyhedraHistorical* pars = GetShape()->GetOriginalParameters();
  G4PolyhedraHistorical* pdata = new G4PolyhedraHistorical(pars->fNumZPlanes);
  pdata->Start_angle = pars->fStartAngle;
  pdata->Opening_angle = pars->fOpeningAngle;
  pdata->numSide = pars->fNumSide;
  for (G4int i=0; i<pars->fNumZPlanes; ++i)
  {
    pdata->Z_values[i] = pars->fZValues[i];
    pdata->Rmin[i] = pars->Rmin[i];
    pdata->Rmax[i] = pars->Rmax[i];
  }
  return pdata;
}
inline void G4UPolyhedra::SetOriginalParameters(G4PolyhedraHistorical* pars)
{
  UPolyhedraHistorical* pdata = GetShape()->GetOriginalParameters();
  pdata->fStartAngle = pars->Start_angle;
  pdata->fOpeningAngle = pars->Opening_angle;
  pdata->fNumSide = pars->numSide;
  pdata->fNumZPlanes = pars->Num_z_planes;
  for (G4int i=0; i<pdata->fNumZPlanes; ++i)
  {
    pdata->fZValues[i] = pars->Z_values[i];
    pdata->Rmin[i] = pars->Rmin[i];
    pdata->Rmax[i] = pars->Rmax[i];
  }
  fRebuildPolyhedron = true;
}
inline G4bool G4UPolyhedra::Reset()
{
  return GetShape()->Reset();
}

#endif
