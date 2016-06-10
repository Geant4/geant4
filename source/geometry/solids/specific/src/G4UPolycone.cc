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
// Implementation of G4UPolycone wrapper class
// --------------------------------------------------------------------

#include "G4Polycone.hh"
#include "G4UPolycone.hh"
#include "G4VPVParameterisation.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor (GEANT3 style parameters)
//  
G4UPolycone::G4UPolycone( const G4String& name, 
                              G4double phiStart,
                              G4double phiTotal,
                              G4int numZPlanes,
                        const G4double zPlane[],
                        const G4double rInner[],
                        const G4double rOuter[]  )
  : G4USolid(name,  new UPolycone(name, phiStart, phiTotal,
                                  numZPlanes, zPlane, rInner, rOuter))
{
}


////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UPolycone::G4UPolycone(const G4String& name, 
                               G4double phiStart,
                               G4double phiTotal,
                               G4int    numRZ,
                         const G4double r[],
                         const G4double z[]   )
  : G4USolid(name, new UPolycone(name, phiStart, phiTotal, numRZ, r, z))
{ 
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UPolycone::G4UPolycone( __void__& a )
  : G4USolid(a)
{
}


////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UPolycone::~G4UPolycone()
{
}


////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UPolycone::G4UPolycone( const G4UPolycone &source )
  : G4USolid( source )
{
}


////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UPolycone &G4UPolycone::operator=( const G4UPolycone &source )
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}


////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4UPolycone::ComputeDimensions(G4VPVParameterisation* p,
                                    const G4int n,
                                    const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Polycone*)this,n,pRep);
}


//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UPolycone::Clone() const
{
  return new G4UPolycone(*this);
}


////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UPolycone::CreatePolyhedron() const
{
  G4PolyconeHistorical* original_parameters = GetOriginalParameters();
  G4PolyhedronPcon*
  polyhedron = new G4PolyhedronPcon( original_parameters->Start_angle,
                                     original_parameters->Opening_angle,
                                     original_parameters->Num_z_planes,
                                     original_parameters->Z_values,
                                     original_parameters->Rmin,
                                     original_parameters->Rmax );

  delete original_parameters;  // delete local copy 

  return polyhedron;
}
