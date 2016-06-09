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
// $Id: G4GeomTestErrorList.cc,v 1.2 2003/11/03 17:15:21 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestErrorList
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4GeomTestErrorList.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHistoryHandle.hh"

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
// Calculate rotation & translation to global coordinates
//
void G4GeomTestErrorList::FindGlobalCoordinateSystem()
{
  G4TouchableHistoryHandle aTouchable =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking()->CreateTouchableHistoryHandle();
  G4AffineTransform globTransform =
    aTouchable->GetHistory()->GetTopTransform().Inverse();

  globalTranslation = globTransform.NetTranslation();
  globalRotation = globTransform.NetRotation();
}
