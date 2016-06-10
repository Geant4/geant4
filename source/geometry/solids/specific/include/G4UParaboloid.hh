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
// G4UParaboloid
//
// Class description:
//
//   Wrapper class for UParaboloid to make use of it from USolids module.

// History:
// 19.08.15 Guilherme Lima, FNAL
// --------------------------------------------------------------------
#ifndef G4UPARABOLOID_HH
#define G4UPARABOLOID_HH

#include "G4USolid.hh"

#if defined(G4GEOM_USE_USOLIDS)

#include "UParaboloid.hh"
#include "G4Polyhedron.hh"

class G4UParaboloid : public G4USolid {

public:  // with description

  G4UParaboloid(const G4String& name, G4double dz, G4double rlo, G4double rhi);

  ~G4UParaboloid();

  G4VSolid* Clone() const;

  inline UParaboloid* GetShape() const;

  inline G4double GetZHalfLength() const;
  inline G4double GetRadiusMinusZ() const;
  inline G4double GetRadiusPlusZ() const;

  // inline G4bool Reset();

public:  // without description

  G4UParaboloid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  G4UParaboloid( const G4UParaboloid &source );
  G4UParaboloid &operator=( const G4UParaboloid &source );
      // Copy constructor and assignment operator.
  // G4Polyhedron* CreatePolyhedron() const;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UParaboloid* G4UParaboloid::GetShape() const
{
  return (UParaboloid*) fShape;
}

inline G4double G4UParaboloid::GetZHalfLength() const {
  return GetShape()->GetDz();
}
inline G4double G4UParaboloid::GetRadiusMinusZ() const {
  return GetShape()->GetRlo();
}
inline G4double G4UParaboloid::GetRadiusPlusZ() const {
  return GetShape()->GetRhi();
}
// inline G4ParaboloidHistorical* G4UParaboloid::GetOriginalParameters() const
// {
//   UParaboloidHistorical* pars = GetShape()->GetOriginalParameters();
//   G4ParaboloidHistorical* pdata = new G4ParaboloidHistorical(pars->fNumZPlanes);
//   pdata->Start_angle = pars->fStartAngle;
//   pdata->Opening_angle = pars->fOpeningAngle;
//   for (G4int i=0; i<pars->fNumZPlanes; ++i)
//   {
//     pdata->Z_values[i] = pars->fZValues[i];
//     pdata->Rmin[i] = pars->Rmin[i];
//     pdata->Rmax[i] = pars->Rmax[i];
//   }
//   return pdata;
// }
// inline void G4UParaboloid::SetOriginalParameters(G4ParaboloidHistorical* pars)
// {
//   UParaboloidHistorical* pdata = GetShape()->GetOriginalParameters();
//   pdata->fStartAngle = pars->Start_angle;
//   pdata->fOpeningAngle = pars->Opening_angle;
//   pdata->fNumZPlanes = pars->Num_z_planes;
//   for (G4int i=0; i<pdata->fNumZPlanes; ++i)
//   {
//     pdata->fZValues[i] = pars->Z_values[i];
//     pdata->Rmin[i] = pars->Rmin[i];
//     pdata->Rmax[i] = pars->Rmax[i];
//   }
//   fRebuildPolyhedron = true;
// }
// inline G4bool G4UParaboloid::Reset()
// {
//   GetShape()->Reset();
//   return 0;
// }

#endif  // G4GEOM_USE_USOLIDS

#endif
