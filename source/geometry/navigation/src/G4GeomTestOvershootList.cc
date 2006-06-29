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
// $Id: G4GeomTestOvershootList.cc,v 1.3 2006-06-29 18:36:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestOvershootList
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

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
{}


//
// Default constructor
//
G4GeomTestOvershootList::G4GeomTestOvershootList()
  : G4GeomTestErrorList(0),
    daughter(0)
{}


//
// Destructor
//
G4GeomTestOvershootList::~G4GeomTestOvershootList()
{}


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
