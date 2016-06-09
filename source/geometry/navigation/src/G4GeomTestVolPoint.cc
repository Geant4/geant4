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
// G4GeomTestVolPoint
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4GeomTestVolPoint.hh"

//
// Constructor (specific)
//
G4GeomTestVolPoint::G4GeomTestVolPoint( const G4ThreeVector &thePoint,
                                              G4double theS,
                                              G4bool isEntering,
                                              G4int theDaughterIndex )
  : G4GeomTestPoint( thePoint, theS, isEntering ),
    daughterIndex(theDaughterIndex)
{;}


//
// Constructor (from base)
//
G4GeomTestVolPoint::G4GeomTestVolPoint( const G4GeomTestPoint &base,
                                              G4int theDaughterIndex )
  : G4GeomTestPoint( base ), daughterIndex(theDaughterIndex)
{;}


//
// Constructor (from base, with coordinate transformation)
//
G4GeomTestVolPoint::G4GeomTestVolPoint( const G4GeomTestPoint &base,
                                              G4int theDaughterIndex,
                                        const G4ThreeVector &translation,
                                        const G4RotationMatrix *rotation )
  : G4GeomTestPoint( base ), daughterIndex(theDaughterIndex)
{
  //
  // Rotate point
  //
  if (rotation)
    p = rotation->inverse()*p - translation;
  else
    p = p - translation;
}


//
// Constructor (copy)
//
G4GeomTestVolPoint::G4GeomTestVolPoint( const G4GeomTestVolPoint &other )
  : G4GeomTestPoint( other ), daughterIndex(other.daughterIndex)
{;}


//
// Constructor (default)
//
G4GeomTestVolPoint::G4GeomTestVolPoint()
  : daughterIndex(-1)
{;}


//
// Destructor
//
G4GeomTestVolPoint::~G4GeomTestVolPoint() {;}


//
// Assignment operator
//
G4GeomTestVolPoint&
G4GeomTestVolPoint::operator=(const G4GeomTestVolPoint& other)
{
   // Check assignment to self
   //
   if (this == &other)  { return *this; }

   // Copy base class data
   //
   G4GeomTestPoint::operator=(other);

   // Copy data
   //
   daughterIndex = other.daughterIndex;

   return *this;
}


//
// Volume accessor
//
G4int G4GeomTestVolPoint::GetDaughterIndex() const
{
  return daughterIndex;
}
