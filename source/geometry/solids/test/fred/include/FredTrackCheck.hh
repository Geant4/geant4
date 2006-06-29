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
// FredTrackCheck.hh
//
// Implementation of Fred's track hit checker
//

#ifndef FredTrackCheck_hh
#define FredTrackCheck_hh

#include "globals.hh"
#include "G4Track.hh"

#include <vector>


class FredTrackData {
	public:
	FredTrackData( G4int thisId ) { id = thisId; nIn = 0; nOut = 0; }
	~FredTrackData() {;}
	
	G4bool operator==(const FredTrackData &right) const { return id == right.id; }
        G4bool operator<(const FredTrackData &right) const { return id < right.id; }
	G4bool operator<=(const FredTrackData &right) const { return id <= right.id; }
	
	void IncrementIn() {nIn++;}
	void IncrementOut() {nOut++;}
	G4int GetNIn() { return nIn; }
	G4int GetNOut() { return nOut; }
	
	private:
	G4int	id, nIn, nOut;
};


class FredTrackCheck {
	
	public:
	FredTrackCheck();
	~FredTrackCheck();
	
	void AddInHit( const G4Track *track );
	void AddOutHit( const G4Track *track );
	
	void GetTrackStat( const G4Track *track, G4int *nIn, G4int *nOut );
	void GetTrackStatID( G4int trackID, G4int *nIn, G4int *nOut );
	void Clear( );

	private:
	std::vector<FredTrackData*> hitList;
};

#endif
