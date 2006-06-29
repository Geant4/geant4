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
// $Id: G4StackedTrack.cc,v 1.9 2006-06-29 18:10:21 gunter Exp $
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

const G4StackedTrack & G4StackedTrack::operator=(const G4StackedTrack &)
{ return *this; }
G4int G4StackedTrack::operator==(const G4StackedTrack &right) const 
{ return (this==&right); }
G4int G4StackedTrack::operator!=(const G4StackedTrack &right) const 
{ return (this!=&right); }


