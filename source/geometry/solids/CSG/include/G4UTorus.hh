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
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4UTorus
//
// Class description:
//
// Wrapper class for UTorus to make use of UTorus from USolids module.

// History:
// 19.08.15 Guilherme Lima, FNAL
// --------------------------------------------------------------------
#ifndef G4UTORUS_HH
#define G4UTORUS_HH

#include "G4USolid.hh"

#if defined(G4GEOM_USE_USOLIDS)

#include "UTorus.hh"
#include "G4Polyhedron.hh"

class G4UTorus : public G4USolid
{
  public:  // with description

    G4UTorus(const G4String& pName,
                   G4double rmin, G4double rmax, G4double rtor,
                   G4double sphi, G4double dphi);
      // Constructs a torus with name and geometrical parameters

   ~G4UTorus();

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline UTorus* GetShape() const;

    inline G4double GetRmin() const;
    inline G4double GetRmax() const;
    inline G4double GetRtor() const;
    inline G4double GetSPhi() const;
    inline G4double GetDPhi() const;

    inline void SetRmin(G4double arg);
    inline void SetRmax(G4double arg);
    inline void SetRtor(G4double arg);
    inline void SetSPhi(G4double arg);
    inline void SetDPhi(G4double arg);

    inline void SetAllParameters(G4double arg1, G4double arg2,
                                 G4double arg3, G4double arg4, G4double arg5);

  public:  // without description

    G4UTorus(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UTorus(const G4UTorus& rhs);
    G4UTorus& operator=(const G4UTorus& rhs);
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UTorus* G4UTorus::GetShape() const
{
  return (UTorus*) fShape;
}

inline G4double G4UTorus::GetRmin() const
{
  return GetShape()->GetRmin();
}

inline G4double G4UTorus::GetRmax() const
{
  return GetShape()->GetRmax();
}

inline G4double G4UTorus::GetRtor() const
{
  return GetShape()->GetRtor();
}

inline G4double G4UTorus::GetSPhi() const
{
  return GetShape()->GetSPhi();
}

inline G4double G4UTorus::GetDPhi() const
{
  return GetShape()->GetDPhi();
}

inline void G4UTorus::SetRmin(G4double arg)
{
  GetShape()->SetRmin(arg);
  fRebuildPolyhedron = true;
}

inline void G4UTorus::SetRmax(G4double arg)
{
  GetShape()->SetRmax(arg);
  fRebuildPolyhedron = true;
}

inline void G4UTorus::SetRtor(G4double arg)
{
  GetShape()->SetRtor(arg);
  fRebuildPolyhedron = true;
}

inline void G4UTorus::SetSPhi(G4double arg)
{
  GetShape()->SetSPhi(arg);
  fRebuildPolyhedron = true;
}

inline void G4UTorus::SetDPhi(G4double arg)
{
  GetShape()->SetDPhi(arg);
  fRebuildPolyhedron = true;
}

inline void G4UTorus::SetAllParameters(G4double arg1, G4double arg2,
                        G4double arg3, G4double arg4, G4double arg5)
{
  GetShape()->SetRmin(arg1);
  GetShape()->SetRmax(arg2);
  GetShape()->SetRtor(arg3);
  GetShape()->SetSPhi(arg4);
  GetShape()->SetDPhi(arg5);
  fRebuildPolyhedron = true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
