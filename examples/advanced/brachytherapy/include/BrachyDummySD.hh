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
//    ********************************
//    *                              *
//    *       BrachyDummySD.hh       *
//    *                              *
//    ********************************

// Dummy sensitive used only to flag sensitivity in cells of RO geometry.

#ifndef BrachyDummySD_h
#define BrachyDummySD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;

class BrachyDummySD : public G4VSensitiveDetector
{
 public:
 	BrachyDummySD();
  	~BrachyDummySD() {};
  
	void Initialize(G4HCofThisEvent*HCE) {};
	G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {return false;}
	void EndOfEvent(G4HCofThisEvent*HCE) {};
	void clear() {};
	void DrawAll() {};
	void PrintAll() {};
};

BrachyDummySD::BrachyDummySD() : G4VSensitiveDetector("dummySD")
{}
#endif
