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
#ifndef G4_CASCADE_SIGMAZEROP_CHANNEL_HH
#define G4_CASCADE_SIGMAZEROP_CHANNEL_HH

#include "G4CascadeChannel.hh"


class G4CascadeSigmaZeroPChannel : public G4CascadeChannel {

public:

  G4CascadeSigmaZeroPChannel();
  virtual ~G4CascadeSigmaZeroPChannel();

  G4double getCrossSection(G4double ke) const; 
  G4int getMultiplicity(G4double ke) const;
  std::vector<G4int> getOutgoingParticleTypes(G4int mult, G4double ke) const;

private:

  static G4double s0ptot[31];
  static G4double s0pMultiplicities[6][31];

  static const G4int s0pindex[6][2];
  static const G4int s0p2bfs[3][2];
  static const G4int s0p3bfs[12][3];
  static const G4int s0p4bfs[33][4];
  static const G4int s0p5bfs[59][5];
  static const G4int s0p6bfs[30][6];
  static const G4int s0p7bfs[20][7];

  static const G4float s0pCrossSections[157][31];
};        

#endif
