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
// Implementation for G4UCons wrapper class
// --------------------------------------------------------------------

#include "G4Cons.hh"
#include "G4UCons.hh"
#include "G4VPVParameterisation.hh"
 
//////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//               - note if pDPhi>2PI then reset to 2PI

G4UCons::G4UCons( const G4String& pName,
                        G4double  pRmin1, G4double pRmax1,
                        G4double  pRmin2, G4double pRmax2,
                        G4double pDz,
                        G4double pSPhi, G4double pDPhi)
  : G4USolid(pName, new UCons(pName, pRmin1, pRmax1, pRmin2, pRmax2,
                                     pDz, pSPhi, pDPhi))
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UCons::G4UCons( __void__& a )
  : G4USolid(a)
{
}

///////////////////////////////////////////////////////////////////////
//
// Destructor

G4UCons::~G4UCons()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UCons::G4UCons(const G4UCons& rhs)
  : G4USolid(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UCons& G4UCons::operator = (const G4UCons& rhs) 
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
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UCons::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int                  n,
                                const G4VPhysicalVolume*     pRep    )
{
  p->ComputeDimensions(*(G4Cons*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UCons::Clone() const
{
  return new G4UCons(*this);
}
