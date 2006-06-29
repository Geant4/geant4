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
#ifndef G4_CASCADE_SIGMAMINUSN_CHANNEL_HH
#define G4_CASCADE_SIGMAMINUSN_CHANNEL_HH

#include "G4CascadeChannel.hh"


class G4CascadeSigmaMinusNChannel : public G4CascadeChannel {

public:

  G4CascadeSigmaMinusNChannel();
  virtual ~G4CascadeSigmaMinusNChannel();

  G4double getCrossSection(G4double ke) const; 
  G4int getMultiplicity(G4double ke) const;
  std::vector<G4int> getOutgoingParticleTypes(G4int mult, G4double ke) const;

private:

  static G4double smntot[31];
  static G4double smnMultiplicities[6][31];

  static const G4int smnindex[6][2];
  static const G4int smn2bfs[1][2];
  static const G4int smn3bfs[6][3];
  static const G4int smn4bfs[20][4];
  static const G4int smn5bfs[42][5];
  static const G4int smn6bfs[25][6];
  static const G4int smn7bfs[17][7];

  static const G4float smnCrossSections[111][31];
};        

#endif
