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
// $Id: G4VRangeTest.hh,v 1.2 2003-01-22 18:42:23 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 5  Oct 2001   MGP        Created
// 22 Jan 2003 V.Ivanchenko Cut per region
//
// -------------------------------------------------------------------
// Class description:
// Base class for a strategy pattern encapsulating algorithms to test the range
// of a particle
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4VRANGETEST_HH
#define G4VRANGETEST_HH 1

#include "globals.hh"
class G4ParticleDefinition;
class G4MaterialCutsCouple;

class G4VRangeTest {

public:

  G4VRangeTest() { }

  virtual ~G4VRangeTest() { }

  virtual G4bool Escape(const G4ParticleDefinition* particle,
			const G4MaterialCutsCouple* couple,
			G4double energy,
			G4double safety) const = 0;

private:

};

#endif











