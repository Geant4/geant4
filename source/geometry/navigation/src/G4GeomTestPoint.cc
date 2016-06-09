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
// $Id: G4GeomTestPoint.cc,v 1.2 2003/11/03 17:15:21 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
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

