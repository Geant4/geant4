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
#ifndef G4_CASCADE_XIMINUSN_CHANNEL_HH
#define G4_CASCADE_XIMINUSN_CHANNEL_HH

#include "G4CascadeChannel.hh"


class G4CascadeXiMinusNChannel : public G4CascadeChannel {

public:

  G4CascadeXiMinusNChannel();
  virtual ~G4CascadeXiMinusNChannel();

  G4double getCrossSection(G4double ke) const; 
  G4int getMultiplicity(G4double ke) const;
  std::vector<G4int> getOutgoingParticleTypes(G4int mult, G4double ke) const;

private:

  static G4double xmntot[31];
  static G4double xmnMultiplicities[6][31];

  static const G4int xmnindex[6][2];
  static const G4int xmn2bfs[3][2];
  static const G4int xmn3bfs[18][3];
  static const G4int xmn4bfs[53][4];
  static const G4int xmn5bfs[2][5];
  static const G4int xmn6bfs[2][6];
  static const G4int xmn7bfs[2][7];

  static const G4float xmnCrossSections[80][31];
};        

#endif
