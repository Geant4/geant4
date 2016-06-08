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
// $Id: G4GeomTestErrorList.cc,v 1.1 2001/10/17 12:59:57 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestErrorList
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)

#include "G4GeomTestErrorList.hh"
#include "G4VPhysicalVolume.hh"

//
// Constructor
//
G4GeomTestErrorList::G4GeomTestErrorList( const G4VPhysicalVolume *theMother )
  : mother(theMother)
{
  FindGlobalCoordinateSystem();
}


//
// Destructor
//
G4GeomTestErrorList::~G4GeomTestErrorList() {;}


//
// AddError
//
// Add a new Error point to our list, where the point
// is given in the coordinate system of the mother
//
void G4GeomTestErrorList::AddError( const G4ThreeVector &s1,
                                    const G4ThreeVector &s2  )
{
  segments.push_back( Segment(s1,s2) );
}


//
// Accessors
//
const G4VPhysicalVolume *G4GeomTestErrorList::GetMother() const
{
  return mother;
}


//
// Return number of segments
//
G4int G4GeomTestErrorList::NumError() const
{
  return segments.size();
}


//
// GetMotherPoints
//
// Return start and end points in the coordinate system
// of the mother
//
void G4GeomTestErrorList::GetMotherPoints( G4int i, 
                                           G4ThreeVector &s1, 
                                           G4ThreeVector &s2 ) const
{
  s1 = segments[i].GetS1();
  s2 = segments[i].GetS2();
}


//
// GetGlobalPoints
//
// Return start and end points in the global coordinate system
//
void G4GeomTestErrorList::GetGlobalPoints( G4int i, 
                                           G4ThreeVector &s1, 
                                           G4ThreeVector &s2 ) const
{
  s1 = globalTranslation + globalRotation*segments[i].GetS1();
  s2 = globalTranslation + globalRotation*segments[i].GetS2();
}


//
// GetOneDaughtPoints
//
// Return start and end points in the coordinate system of
// a daughter
//
void G4GeomTestErrorList::GetOneDaughtPoints( const G4VPhysicalVolume *daught,
                                                    G4int i, 
                                                    G4ThreeVector &s1, 
                                                    G4ThreeVector &s2 ) const
{
  const G4RotationMatrix *rotation = daught->GetFrameRotation();
  const G4ThreeVector &translation = daught->GetFrameTranslation();

  if (rotation) {
    s1 = (*rotation)*(translation + segments[i].GetS1());
    s2 = (*rotation)*(translation + segments[i].GetS2());
  }
  else {
    s1 = translation + segments[i].GetS1();
    s2 = translation + segments[i].GetS2();
  }
}



//
// FindGlobalCoordinateSystem
//
// Calculate translation to global coordinates
//
void G4GeomTestErrorList::FindGlobalCoordinateSystem()
{
  globalTranslation = G4ThreeVector(0,0,0);
  globalRotation = G4RotationMatrix();
  
  const G4VPhysicalVolume *vol = mother;
  
  for(;;) {
    const G4VPhysicalVolume *parent = vol->GetMother();
    if (parent == 0) break;

    const G4RotationMatrix &rotation = vol->GetObjectRotationValue();
    const G4ThreeVector &translation = vol->GetObjectTranslation();
    
    globalTranslation = translation + rotation*globalTranslation;
    globalRotation = rotation*globalRotation;
  
    vol = parent;
  }
}
