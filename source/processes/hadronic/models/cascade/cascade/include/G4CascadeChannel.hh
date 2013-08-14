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
// $Id: G4CascadeChannel.hh 67796 2013-03-08 06:18:39Z mkelsey $
//
// 20100514  M. Kelsey -- All functionality removed except quantum-number
//		validation functions.
// 20110719  M. Kelsey -- Repurpose as abstract base class for all channels.
// 20110802  M. Kelsey -- Fix Coverity #29838: Virtual dtor required for base.
// 20110923  M. Kelsey -- Add optional ostream& argument to printTable(), add
//		stream operator<<, implemented in new .cc file.

#ifndef G4_CASCADE_CHANNEL_HH
#define G4_CASCADE_CHANNEL_HH

#include "globals.hh"
#include <iosfwd>
#include <vector>


class G4CascadeChannel {
public:
  G4CascadeChannel() {}
  virtual ~G4CascadeChannel() {}

  virtual G4double getCrossSection(double ke) const = 0;
  virtual G4double getCrossSectionSum(double ke) const = 0;
  virtual G4int getMultiplicity(G4double ke) const = 0;

  virtual void getOutgoingParticleTypes(std::vector<G4int>& kinds,
					G4int mult, G4double ke) const = 0;

  virtual void printTable(std::ostream& os=G4cout) const = 0;
};

std::ostream& operator<<(std::ostream& os, const G4CascadeChannel& chan);

#endif	/* G4_CASCADE_CHANNEL_HH */
