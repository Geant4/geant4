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
#ifndef G4_CASCADE_KPLUSP_CHANNEL_HH
#define G4_CASCADE_KPLUSP_CHANNEL_HH

#include "G4CascadeChannel.hh"


class G4CascadeKplusPChannel : public G4CascadeChannel {

public:

  G4CascadeKplusPChannel();
  virtual ~G4CascadeKplusPChannel();

  G4double getCrossSection(G4double ke) const; 
  G4int getMultiplicity(G4double ke) const;
  std::vector<G4int> getOutgoingParticleTypes(G4int mult, G4double ke) const;

private:

  static G4double kpptot[31];
  static G4double kppMultiplicities[6][31];

  static const G4int kppindex[6][2];
  static const G4int kpp2bfs[1][2];
  static const G4int kpp3bfs[4][3];
  static const G4int kpp4bfs[10][4];
  static const G4int kpp5bfs[19][5];
  static const G4int kpp6bfs[28][6];
  static const G4int kpp7bfs[38][7];

  static const G4float kppCrossSections[100][31];
};        

#endif
