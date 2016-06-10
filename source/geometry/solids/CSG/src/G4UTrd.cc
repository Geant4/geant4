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
// Implementation for G4UTrd wrapper class
// --------------------------------------------------------------------

#include "G4Trd.hh"
#include "G4UTrd.hh"
#include "G4VPVParameterisation.hh"

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths
//
G4UTrd::G4UTrd(const G4String& pName,
                     G4double pdx1,  G4double pdx2,
                     G4double pdy1,  G4double pdy2,
                     G4double pdz)
  : G4USolid(pName, new UTrd(pName, pdx1, pdx2, pdy1, pdy2, pdz))
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTrd::G4UTrd( __void__& a )
  : G4USolid(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UTrd::~G4UTrd()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTrd::G4UTrd(const G4UTrd& rhs)
  : G4USolid(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTrd& G4UTrd::operator = (const G4UTrd& rhs) 
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
//
void G4UTrd::ComputeDimensions(      G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Trd*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UTrd::Clone() const
{
  return new G4UTrd(*this);
}
