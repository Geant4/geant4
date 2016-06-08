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
// $Id: G4He3CoulombBarrier.hh,v 1.1.8.1 2001/06/28 19:13:02 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4He3CoulombBarrier_h
#define G4He3CoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4He3CoulombBarrier : public G4CoulombBarrier
{
public:
	G4He3CoulombBarrier() : G4CoulombBarrier(3,2) {};
	~G4He3CoulombBarrier() {};

private:
	G4He3CoulombBarrier(const G4He3CoulombBarrier & right);

	const G4He3CoulombBarrier & operator=(const G4He3CoulombBarrier & right);
	G4bool operator==(const G4He3CoulombBarrier & right) const;
	G4bool operator!=(const G4He3CoulombBarrier & right) const;
  
private:

	virtual G4double BarrierPenetrationFactor(const G4double aZ) const;


};

#endif
