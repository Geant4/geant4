//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4GeomTestOvershootList.cc,v 1.1 2001/10/17 12:59:58 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestOvershootList
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)

#include "G4GeomTestOvershootList.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

//
// Constructor
//
G4GeomTestOvershootList::G4GeomTestOvershootList(
                             const G4VPhysicalVolume *theMother,
                                   G4int theDaughterIndex  )
  : G4GeomTestErrorList(theMother),
    daughter(theDaughterIndex)
{;}


//
// Default constructor
//
G4GeomTestOvershootList::G4GeomTestOvershootList()
  : G4GeomTestErrorList(0),
    daughter(0)
{;}


//
// Destructor
//
G4GeomTestOvershootList::~G4GeomTestOvershootList() {;}


//
// Comparison operators
//
G4bool
G4GeomTestOvershootList::operator==( const G4GeomTestOvershootList &other ) const
{
  return daughter==other.daughter;
}

G4bool
G4GeomTestOvershootList::operator< ( const G4GeomTestOvershootList &other ) const
{
  return (daughter < other.daughter);
}



//
// Accessors
//
const G4VPhysicalVolume *G4GeomTestOvershootList::GetDaughter() const
{ 
  return GetMother()->GetLogicalVolume()->GetDaughter(daughter);
}

G4int G4GeomTestOvershootList::GetDaughterIndex() const 
{
  return daughter;
}


//
// GetDaughtPoints
//
// Return start and end points in the coordinate system of
// the daughter
//
void G4GeomTestOvershootList::GetDaughtPoints( G4int i, 
                                               G4ThreeVector &s1, 
                                               G4ThreeVector &s2 ) const
{
  GetOneDaughtPoints( GetDaughter(), i, s1, s2 );
}
