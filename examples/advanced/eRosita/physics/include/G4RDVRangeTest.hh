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

#ifndef G4RDVRANGETEST_HH
#define G4RDVRANGETEST_HH 1

#include "globals.hh"
class G4ParticleDefinition;
class G4MaterialCutsCouple;

class G4RDVRangeTest {

public:

  G4RDVRangeTest() { }

  virtual ~G4RDVRangeTest() { }

  virtual G4bool Escape(const G4ParticleDefinition* particle,
			const G4MaterialCutsCouple* couple,
			G4double energy,
			G4double safety) const = 0;

private:

};

#endif











