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
// G4GeomTestPoint
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4GeomTestPoint.hh"


//
// Default constructor
//
G4GeomTestPoint::G4GeomTestPoint()
  : p(0),
    s(0),
    entering(false)
{;}


//
// Specific constructor
//
G4GeomTestPoint::G4GeomTestPoint( const G4ThreeVector &thePoint,
                                        G4double theS,
                                        G4bool isEntering )
  : p(thePoint),
    s(theS),
    entering(isEntering)
{;}


//
// Copy constructor
//
G4GeomTestPoint::G4GeomTestPoint( const G4GeomTestPoint &other )
  : p(other.p),
    s(other.s),
    entering(other.entering)
{;}


//
// Destructor
//
G4GeomTestPoint::~G4GeomTestPoint() {;}

//
// Assignment operator
//
G4GeomTestPoint& G4GeomTestPoint::operator=(const G4GeomTestPoint& other)
{
   // Check assignment to self
   //
   if (this == &other)  { return *this; }

   // Copy data
   //
   p = other.p;
   s = other.s;
   entering = other.entering;

   return *this;
}


//
// Equivalence operator
//
G4bool G4GeomTestPoint::operator==( const G4GeomTestPoint &other ) const
{
  return s == other.s;
}


//
// Order operators
//
G4bool G4GeomTestPoint::operator<( const G4GeomTestPoint &other ) const
{
  return s < other.s;
}

G4bool G4GeomTestPoint::operator<=( const G4GeomTestPoint &other ) const
{
  return s <= other.s;
}


//
// Return position
//
const G4ThreeVector &G4GeomTestPoint::GetPosition() const
{
  return p;
}


//
// Return distance
//
G4double G4GeomTestPoint::GetDistance() const
{
  return s;
}


//
// Return true if point was entering
//
G4bool G4GeomTestPoint::Entering() const
{
  return entering;
}

