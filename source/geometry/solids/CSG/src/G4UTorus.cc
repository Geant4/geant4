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
// Implementation for G4UTorus wrapper class
//
// 19-08-2015 Guilherme Lima, FNAL
//
// --------------------------------------------------------------------

#include "G4Torus.hh"
#include "G4UTorus.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4VPVParameterisation.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths


G4UTorus::G4UTorus(const G4String& pName,
                         G4double rmin, G4double rmax, G4double rtor,
                         G4double sphi, G4double dphi)
  : G4USolid(pName, new UTorus(pName, rmin, rmax, rtor, sphi, dphi))
{ }

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UTorus::G4UTorus( __void__& a )
  : G4USolid(a)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UTorus::~G4UTorus() { }

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UTorus::G4UTorus(const G4UTorus& rhs)
  : G4USolid(rhs)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UTorus& G4UTorus::operator = (const G4UTorus& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4USolid::operator=(rhs);

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers

G4double G4UTorus::GetRmin() const
{
  return GetShape()->GetRmin();
}

G4double G4UTorus::GetRmax() const
{
  return GetShape()->GetRmax();
}

G4double G4UTorus::GetRtor() const
{
  return GetShape()->GetRtor();
}

G4double G4UTorus::GetSPhi() const
{
  return GetShape()->GetSPhi();
}

G4double G4UTorus::GetDPhi() const
{
  return GetShape()->GetDPhi();
}

void G4UTorus::SetRmin(G4double arg)
{
  GetShape()->SetRmin(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetRmax(G4double arg)
{
  GetShape()->SetRmax(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetRtor(G4double arg)
{
  GetShape()->SetRtor(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetSPhi(G4double arg)
{
  GetShape()->SetSPhi(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetDPhi(G4double arg)
{
  GetShape()->SetDPhi(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetAllParameters(G4double arg1, G4double arg2,
                        G4double arg3, G4double arg4, G4double arg5)
{
  GetShape()->SetRmin(arg1);
  GetShape()->SetRmax(arg2);
  GetShape()->SetRtor(arg3);
  GetShape()->SetSPhi(arg4);
  GetShape()->SetDPhi(arg5);
  fRebuildPolyhedron = true;
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UTorus::ComputeDimensions(G4VPVParameterisation* p,
                                 const G4int n,
                                 const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Torus*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UTorus::Clone() const
{
  return new G4UTorus(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization

G4Polyhedron* G4UTorus::CreatePolyhedron() const
{
  return new G4PolyhedronTorus(GetRmin(),
                               GetRmax(),
                               GetRtor(),
                               GetSPhi(),
                               GetDPhi());
}

#endif  // G4GEOM_USE_USOLIDS
