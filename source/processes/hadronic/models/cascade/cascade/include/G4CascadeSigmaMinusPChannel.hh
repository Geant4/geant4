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
#ifndef G4_CASCADE_SIGMAMINUSP_CHANNEL_HH
#define G4_CASCADE_SIGMAMINUSP_CHANNEL_HH

#include "G4CascadeChannel.hh"


class G4CascadeSigmaMinusPChannel : public G4CascadeChannel {

public:

  G4CascadeSigmaMinusPChannel();
  virtual ~G4CascadeSigmaMinusPChannel();

  G4double getCrossSection(G4double ke) const; 
  G4int getMultiplicity(G4double ke) const;
  std::vector<G4int> getOutgoingParticleTypes(G4int mult, G4double ke) const;

private:

  static G4double smptot[31];
  static G4double smpMultiplicities[6][31];

  static const G4int smpindex[6][2];
  static const G4int smp2bfs[3][2];
  static const G4int smp3bfs[12][3];
  static const G4int smp4bfs[33][4];
  static const G4int smp5bfs[59][5];
  static const G4int smp6bfs[30][6];
  static const G4int smp7bfs[20][7];

  static const G4float smpCrossSections[157][31];
};        

#endif
