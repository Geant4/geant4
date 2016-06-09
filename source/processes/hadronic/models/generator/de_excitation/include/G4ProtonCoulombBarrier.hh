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
// $Id: G4ProtonCoulombBarrier.hh,v 1.4 2002/12/12 19:17:10 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4ProtonCoulombBarrier_h
#define G4ProtonCoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4ProtonCoulombBarrier : public G4CoulombBarrier
{
public:
	G4ProtonCoulombBarrier() : G4CoulombBarrier(1,1) {};
	~G4ProtonCoulombBarrier() {};

private:
	G4ProtonCoulombBarrier(const G4ProtonCoulombBarrier & right);

	const G4ProtonCoulombBarrier & operator=(const G4ProtonCoulombBarrier & right);
	G4bool operator==(const G4ProtonCoulombBarrier & right) const;
	G4bool operator!=(const G4ProtonCoulombBarrier & right) const;
  
private:

	virtual G4double BarrierPenetrationFactor(const G4double aZ) const;


};

#endif
