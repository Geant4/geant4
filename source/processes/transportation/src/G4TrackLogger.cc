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
// $Id: G4TrackLogger.cc,v 1.1 2002-07-10 15:51:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4TrackLogger.cc
//
// ----------------------------------------------------------------------

#include "G4TrackLogger.hh"

G4TrackLogger::G4TrackLogger() : 
  fPreviousEventID(-1)
{}

G4TrackLogger::~G4TrackLogger(){}


void G4TrackLogger::SetEventID(G4int id){
  if (id != fPreviousEventID) {
    fTrackIDsSet.clear();
    fPreviousEventID =id;
  }
};

G4bool G4TrackLogger::FirstEnterance(G4int trid){
  G4bool first = true;
  G4int n = fTrackIDsSet.count(trid);
  if (n==1) {
    first=false;
    fTrackIDsSet.insert(trid);
  }
  else if (n>1) {
    G4cout << "Error G4TrackLogger::FirstEnterance: " 
	   << "more than one elm in set!" << G4endl;
  }
  return first;
};

