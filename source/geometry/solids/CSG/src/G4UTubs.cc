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
// Implementation for G4UTubs wrapper class
// --------------------------------------------------------------------

#include "G4Tubs.hh"
#include "G4UTubs.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4VPVParameterisation.hh"

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4UTubs::G4UTubs( const G4String& pName,
                        G4double pRMin, G4double pRMax,
                        G4double pDz,
                        G4double pSPhi, G4double pDPhi )
  : G4USolid(pName, new UTubs(pName, pRMin, pRMax, pDz, pSPhi, pDPhi))
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTubs::G4UTubs( __void__& a )
  : G4USolid(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UTubs::~G4UTubs()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UTubs::G4UTubs(const G4UTubs& rhs)
  : G4USolid(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UTubs& G4UTubs::operator = (const G4UTubs& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4USolid::operator=(rhs);

   return *this;
}

/////////////////////////////////////////////////////////////////////////
//
// Accessors and modifiers

G4double G4UTubs::GetInnerRadius() const
{
  return GetShape()->GetInnerRadius();
}
G4double G4UTubs::GetOuterRadius() const
{
  return GetShape()->GetOuterRadius();
}
G4double G4UTubs::GetZHalfLength() const
{
  return GetShape()->GetZHalfLength();
}
G4double G4UTubs::GetStartPhiAngle() const
{
  return GetShape()->GetStartPhiAngle();
}
G4double G4UTubs::GetDeltaPhiAngle() const
{
  return GetShape()->GetDeltaPhiAngle();
}

void G4UTubs::SetInnerRadius(G4double newRMin)
{
  GetShape()->SetInnerRadius(newRMin);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetOuterRadius(G4double newRMax)
{
  GetShape()->SetOuterRadius(newRMax);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetZHalfLength(G4double newDz)
{
  GetShape()->SetZHalfLength(newDz);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetStartPhiAngle(G4double newSPhi, G4bool trig)
{
  GetShape()->SetStartPhiAngle(newSPhi, trig);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetDeltaPhiAngle(G4double newDPhi)
{
  GetShape()->SetDeltaPhiAngle(newDPhi);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UTubs::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*(G4Tubs*)this,n,pRep) ;
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UTubs::Clone() const
{
  return new G4UTubs(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization
//
G4Polyhedron* G4UTubs::CreatePolyhedron() const
{
  return new G4PolyhedronTubs(GetInnerRadius(),
                              GetOuterRadius(),
                              GetZHalfLength(),
                              GetStartPhiAngle(),
                              GetDeltaPhiAngle());
}

#endif  // G4GEOM_USE_USOLIDS
