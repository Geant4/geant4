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
#ifndef Tst34EventAction_h
#define Tst34EventAction_h

#include "globals.hh"
#include "G4UserEventAction.hh"


class Tst34EventAction: public G4UserEventAction {
	public:
	Tst34EventAction();
	~Tst34EventAction();
	
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);
	private:
	G4int nevent;
	G4double dtime;
	G4int calorimeterCollectionId;
};
#endif
