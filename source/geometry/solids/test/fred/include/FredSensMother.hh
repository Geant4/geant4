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
// FredSensMother.hh
//
// Definition of Fred's sensitive mother volume
//

#ifndef FredSensMother_hh
#define FredSensMother_hh

#include "G4VSensitiveDetector.hh"
#include "FredHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class FredSensMother : public G4VSensitiveDetector
{
	public:
	FredSensMother( G4String name );
	~FredSensMother();
	
	//
	// G4V virtual functions all of which must be overridden
	//
	void Initialize( G4HCofThisEvent *HCE );
	G4bool ProcessHits( G4Step *step, G4TouchableHistory *hist );

	//
	// We hold onto the hits here
	//
	private:
	FredHitCollection *hits;
};

#endif
