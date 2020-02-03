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
// Implementation for G4UOrb wrapper class
//
// 30.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------

#include "G4Orb.hh"
#include "G4UOrb.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4VPVParameterisation.hh"
#include "G4PhysicalConstants.hh"

////////////////////////////////////////////////////////////////////////
//
// constructor - check positive radius
//             

G4UOrb::G4UOrb( const G4String& pName, G4double pRmax )
  : G4USolid(pName, new UOrb(pName, pRmax))
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UOrb::G4UOrb( __void__& a )
  : G4USolid(a)
{
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4UOrb::~G4UOrb()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UOrb::G4UOrb(const G4UOrb& rhs)
  : G4USolid(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UOrb& G4UOrb::operator = (const G4UOrb& rhs) 
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

G4double G4UOrb::GetRadius() const
{
  return GetShape()->GetRadius();
}

void G4UOrb::SetRadius(G4double newRmax)
{
  GetShape()->SetRadius(newRmax);
  fRebuildPolyhedron = true;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UOrb::ComputeDimensions(       G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*(G4Orb*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UOrb::Clone() const
{
  return new G4UOrb(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization

G4Polyhedron* G4UOrb::CreatePolyhedron() const
{
  return new G4PolyhedronSphere(0., GetRadius(), 0., twopi, 0., pi);
}

#endif  // G4GEOM_USE_USOLIDS
