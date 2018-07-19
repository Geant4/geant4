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
// $Id: G4RDRangeTest.hh 107396 2017-11-10 08:28:08Z gcosmo $
// GEANT4 tag $Name: geant4-09-01-ref-00 $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 5 Oct  2001   MGP        Created
// 21 Jan 2003   VI         Cuts per region
//
// -------------------------------------------------------------------
// Class description:
// Class for a strategy pattern encapsulating algorithms to test the range
// of a perticle: test range against cut and safety
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4RDRANGETEST_HH
#define G4RDRANGETEST_HH 1

#include "globals.hh"
#include "G4RDVRangeTest.hh"

class G4ParticleDefinition;
class G4Material;

class G4RDRangeTest : public G4RDVRangeTest {

public:

  G4RDRangeTest() { }

  virtual ~G4RDRangeTest();

  virtual G4bool Escape(const G4ParticleDefinition* particle,
			const G4MaterialCutsCouple* couple,
			G4double energy,
			G4double safety) const;

private:

};
 
#endif
 










