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
// $Id$
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestOverlapList
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4GeomTestOverlapList.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

//
// Constructor
//
G4GeomTestOverlapList::G4GeomTestOverlapList(
                             const G4VPhysicalVolume *theMother,
                                   G4int theDaughter1,
                                   G4int theDaughter2 )
  : G4GeomTestErrorList(theMother),
    daughter1(theDaughter1),
    daughter2(theDaughter2)
{}


//
// Default constructor
//
G4GeomTestOverlapList::G4GeomTestOverlapList()
  : G4GeomTestErrorList(0),
    daughter1(0),
    daughter2(0)
{}


//
// Destructor
//
G4GeomTestOverlapList::~G4GeomTestOverlapList()
{}


//
// Comparison operators
//
G4bool
G4GeomTestOverlapList::operator==( const G4GeomTestOverlapList &other ) const
{
  return daughter1==other.daughter1 && daughter2==other.daughter2;
}

G4bool
G4GeomTestOverlapList::operator< ( const G4GeomTestOverlapList &other ) const
{
  if (daughter1 > other.daughter1) return false;
  if (daughter1 < other.daughter1) return true;
  
  return (daughter2 < other.daughter2);
}



//
// Accessors
//
const G4VPhysicalVolume *G4GeomTestOverlapList::GetDaughter1() const
{ 
  return GetMother()->GetLogicalVolume()->GetDaughter(daughter1); 
}

const G4VPhysicalVolume *G4GeomTestOverlapList::GetDaughter2() const
{ 
  return GetMother()->GetLogicalVolume()->GetDaughter(daughter2); 
}

G4int G4GeomTestOverlapList::GetDaughter1Index() const
{
  return daughter1;
}

G4int G4GeomTestOverlapList::GetDaughter2Index() const
{
  return daughter2;
}


//
// GetDaught1Points
//
// Return start and end points in the coordinate system of
// the first daughter
//
void G4GeomTestOverlapList::GetDaught1Points( G4int i, 
                                              G4ThreeVector &s1, 
                                              G4ThreeVector &s2 ) const
{
  GetOneDaughtPoints( GetDaughter1(), i, s1, s2 );
}


//
// GetDaught2Points
//
// Return start and end points in the coordinate system of
// the second daughter
//
void G4GeomTestOverlapList::GetDaught2Points( G4int i, 
                                              G4ThreeVector &s1, 
                                              G4ThreeVector &s2 ) const
{
  GetOneDaughtPoints( GetDaughter2(), i, s1, s2 );
}
