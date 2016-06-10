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
// $Id: G4TrackStack.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
//  Last Modification : 09/Dec/96 M.Asai
//


#ifndef G4TrackStack_h
#define G4TrackStack_h 1

#include "G4StackedTrack.hh"
#include "G4Types.hh"
#include <vector>

class G4SmartTrackStack;

// class description:
//
// This is a stack class used by G4StackManager. This class object
// stores G4StackedTrack class objects in the form of bi-directional
// linked list.

class G4TrackStack : public std::vector<G4StackedTrack>
{
public:
	G4TrackStack() : safetyValve1(0), safetyValve2(0), nstick(0) {}
  G4TrackStack(size_t n) : safetyValve1(4*n/5), safetyValve2(4*n/5-100), nstick(100) { reserve(n);}
  ~G4TrackStack();
  
private:
	const G4TrackStack & operator=(const G4TrackStack &right);
	G4int operator==(const G4TrackStack &right) const;
	G4int operator!=(const G4TrackStack &right) const;
  
public:
	void PushToStack(const G4StackedTrack& aStackedTrack) { push_back(aStackedTrack); }
	G4StackedTrack PopFromStack() { G4StackedTrack st = back(); pop_back(); return st; }
	void TransferTo(G4TrackStack* aStack);
	void TransferTo(G4SmartTrackStack* aStack);
  
        void clearAndDestroy();
private:
	G4int safetyValve1;
  G4int safetyValve2;
	G4int nstick;
  
public:
	G4int GetNTrack() const { return size(); }
	G4int GetMaxNTrack() const { return max_size(); }
  inline G4int GetSafetyValve1() const { return safetyValve1; }
	inline G4int GetSafetyValve2() const { return safetyValve2; }
	inline G4int GetNStick() const { return nstick; }
  
	G4double getTotalEnergy(void) const;
	void SetSafetyValve2(int x) { safetyValve2 = x  < 0 ? 0 : x; }
  
};

#endif
