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
#ifndef G4_CASCADE_KZERON_CHANNEL_HH
#define G4_CASCADE_KZERON_CHANNEL_HH

#include "G4CascadeChannel.hh"


class G4CascadeKzeroNChannel : public G4CascadeChannel {

public:

  G4CascadeKzeroNChannel();
  virtual ~G4CascadeKzeroNChannel();

  G4double getCrossSection(G4double ke) const; 
  G4int getMultiplicity(G4double ke) const;
  std::vector<G4int> getOutgoingParticleTypes(G4int mult, G4double ke) const;

private:

  static G4double k0ntot[31];
  static G4double k0nMultiplicities[6][31];

  static const G4int k0nindex[6][2];
  static const G4int k0n2bfs[1][2];
  static const G4int k0n3bfs[4][3];
  static const G4int k0n4bfs[10][4];
  static const G4int k0n5bfs[19][5];
  static const G4int k0n6bfs[28][6];
  static const G4int k0n7bfs[38][7];

  static const G4float k0nCrossSections[100][31];
};        

#endif
