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
// G4UPolycone
//
// Class description:
//
//   Wrapper class for UPolycone to make use of it from USolids module.

// History:
// 31.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UPOLYCONE_HH
#define G4UPOLYCONE_HH

#include "G4USolid.hh"
#include "UPolycone.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyconeHistorical.hh"
#include "G4Polyhedron.hh"

class G4UPolycone : public G4USolid 
{

  public:  // with description

    G4UPolycone(const G4String& name, 
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int numZPlanes,     // number of z planes
                const G4double zPlane[],    // position of z planes
                const G4double rInner[],    // tangent distance to inner surface
                const G4double rOuter[]  ); // tangent distance to outer surface

    G4UPolycone(const G4String& name, 
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int    numRZ,       // number corners in r,z space
                const G4double r[],         // r coordinate of these corners
                const G4double z[]       ); // z coordinate of these corners

   ~G4UPolycone();
  
    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline UPolycone* GetShape() const;

    inline G4double GetStartPhi()  const;
    inline G4double GetEndPhi()    const;
    inline G4bool IsOpen()         const;
    inline G4int  GetNumRZCorner() const;
    inline G4PolyconeSideRZ GetCorner(G4int index) const;
    inline G4PolyconeHistorical* GetOriginalParameters() const;
    inline void SetOriginalParameters(G4PolyconeHistorical* pars);

    inline G4bool Reset();

  public:  // without description

    G4UPolycone(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UPolycone( const G4UPolycone &source );
    G4UPolycone &operator=( const G4UPolycone &source );
      // Copy constructor and assignment operator.
    G4Polyhedron* CreatePolyhedron() const;
    
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UPolycone* G4UPolycone::GetShape() const
{
  return (UPolycone*) fShape;
}

inline G4double G4UPolycone::GetStartPhi() const
{
  return GetShape()->GetStartPhi();
}
inline G4double G4UPolycone::GetEndPhi() const
{
  return GetShape()->GetEndPhi();
}
inline G4bool G4UPolycone::IsOpen() const
{
  return GetShape()->IsOpen();
}
inline G4int G4UPolycone::GetNumRZCorner() const
{
  return GetShape()->GetNumRZCorner();
}
inline G4PolyconeSideRZ G4UPolycone::GetCorner(G4int index) const
{
  UPolyconeSideRZ pside = GetShape()->GetCorner(index);
  G4PolyconeSideRZ psiderz = { pside.r, pside.z };

  return psiderz;
}
inline G4PolyconeHistorical* G4UPolycone::GetOriginalParameters() const
{
  UPolyconeHistorical* pars = GetShape()->GetOriginalParameters();
  G4PolyconeHistorical* pdata = new G4PolyconeHistorical(pars->fNumZPlanes);
  pdata->Start_angle = pars->fStartAngle;
  pdata->Opening_angle = pars->fOpeningAngle;
  for (G4int i=0; i<pars->fNumZPlanes; ++i)
  {
    pdata->Z_values[i] = pars->fZValues[i];
    pdata->Rmin[i] = pars->Rmin[i];
    pdata->Rmax[i] = pars->Rmax[i];
  }
  return pdata;
}
inline void G4UPolycone::SetOriginalParameters(G4PolyconeHistorical* pars)
{
  UPolyconeHistorical* pdata = GetShape()->GetOriginalParameters();
  pdata->fStartAngle = pars->Start_angle;
  pdata->fOpeningAngle = pars->Opening_angle;
  pdata->fNumZPlanes = pars->Num_z_planes;
  for (G4int i=0; i<pdata->fNumZPlanes; ++i)
  {
    pdata->fZValues[i] = pars->Z_values[i];
    pdata->Rmin[i] = pars->Rmin[i];
    pdata->Rmax[i] = pars->Rmax[i];
  }
  fRebuildPolyhedron = true;
}
inline G4bool G4UPolycone::Reset()
{
  GetShape()->Reset();
  return 0;
}

#endif
