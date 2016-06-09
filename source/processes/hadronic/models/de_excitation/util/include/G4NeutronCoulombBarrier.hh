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
// $Id: G4NeutronCoulombBarrier.hh,v 1.2 2005/06/04 13:29:20 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4NeutronCoulombBarrier_h
#define G4NeutronCoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4NeutronCoulombBarrier : public G4CoulombBarrier
{
public:
	G4NeutronCoulombBarrier() : G4CoulombBarrier(1,0) {};
	~G4NeutronCoulombBarrier() {};

private:
	G4NeutronCoulombBarrier(const G4NeutronCoulombBarrier & right);

	const G4NeutronCoulombBarrier & operator=(const G4NeutronCoulombBarrier & right);
	G4bool operator==(const G4NeutronCoulombBarrier & right) const;
	G4bool operator!=(const G4NeutronCoulombBarrier & right) const;
  
private:

	virtual G4double BarrierPenetrationFactor(const G4double ) const
	{ return 1.0;}


};

#endif
