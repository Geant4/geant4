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
// $Id: G4StackedTrack.cc,v 1.5 2001-07-13 15:01:55 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 02/Feb/96 M.Asai
//

#include "G4StackedTrack.hh"

G4Allocator<G4StackedTrack> aStackedTrackAllocator;

G4StackedTrack::G4StackedTrack() 
:priorityWeight(0.),track(0),trajectory(0),
 previousStackedTrack(0),nextStackedTrack(0)
{ }

G4StackedTrack::G4StackedTrack(G4Track * newTrack, G4VTrajectory * aTrajectory) 
:priorityWeight(0.),track(newTrack),trajectory(aTrajectory),
 previousStackedTrack(0),nextStackedTrack(0)
{ }

G4StackedTrack::~G4StackedTrack()
{ }

const G4StackedTrack & G4StackedTrack::operator=(const G4StackedTrack &right)
{ return *this; }
G4int G4StackedTrack::operator==(const G4StackedTrack &right) const 
{ return false; }
G4int G4StackedTrack::operator!=(const G4StackedTrack &right) const 
{ return true; }


