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
// $Id: G4StackedTrack.cc,v 1.4 2001-07-11 09:58:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 02/Feb/96 M.Asai
//

#include "G4StackedTrack.hh"

G4Allocator<G4StackedTrack> aStackedTrackAllocator;

G4StackedTrack::G4StackedTrack() 
:track(NULL),trajectory(NULL),priorityWeight(0.),
 previousStackedTrack(NULL),nextStackedTrack(NULL)
{ }

G4StackedTrack::G4StackedTrack(G4Track * newTrack, G4VTrajectory * aTrajectory) 
:track(newTrack),trajectory(aTrajectory),priorityWeight(0.),
 previousStackedTrack(NULL),nextStackedTrack(NULL)
{ }

G4StackedTrack::~G4StackedTrack()
{ }

const G4StackedTrack & G4StackedTrack::operator=(const G4StackedTrack &right)
{ return *this; }
int G4StackedTrack::operator==(const G4StackedTrack &right) const 
{ return false; }
int G4StackedTrack::operator!=(const G4StackedTrack &right) const 
{ return true; }


