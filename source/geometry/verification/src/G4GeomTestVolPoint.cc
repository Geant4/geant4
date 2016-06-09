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
// $Id: G4GeomTestVolPoint.cc,v 1.1 2001/10/17 13:00:00 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestVolPoint
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)

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
// Volume accessor
//
G4int G4GeomTestVolPoint::GetDaughterIndex() const
{
  return daughterIndex;
}
