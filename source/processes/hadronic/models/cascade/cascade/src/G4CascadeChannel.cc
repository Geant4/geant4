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

#include "G4CascadeChannel.hh"
#include "Randomize.hh"


G4CascadeChannel::G4CascadeChannel()
{;}
 
G4CascadeChannel::~G4CascadeChannel()
{;}


std::pair<G4int, G4double> 
G4CascadeChannel::interpolateEnergy(G4double e) const
{
  G4int index = 30;
  G4double fraction = 0.0;

  for (G4int i = 1; i < 31; i++) {
    if (e < energyScale[i]) {
      index = i-1;
      fraction = (e - energyScale[index]) / (energyScale[i] - energyScale[index]);
      break;
    }
  }
  return std::pair<G4int, G4double>(index, fraction);
}


G4int 
G4CascadeChannel::sampleFlat(std::vector<G4double> sigma) const
{
  G4int i;
  G4double sum(0.);
  for (i = 0; i < G4int(sigma.size()); i++) sum += sigma[i];
 
  G4double fsum = sum*G4UniformRand();
  G4double partialSum = 0.0;
  G4int channel = 0;

  for (i = 0; i < G4int(sigma.size()); i++) {
    partialSum += sigma[i];
    if (fsum < partialSum) {
      channel = i;
      break;
    }
  }

  return channel;
}


std::vector<G4int> 
G4CascadeChannel::getQnums(G4int type) const
{
  G4int bary=0, str=0, ch=0;
  std::vector<G4int> Qnums(3);

    switch(type) {
    case 1: // proton
      bary = 1;
      str = 0;
      ch = 1;
      break;
    case 2: // neutron
      bary = 1;
      str = 0;
      ch = 0;
      break;
    case 3: // pi+
      bary = 0;
      str = 0;
      ch = 1;
      break;
    case 5: // pi-
      bary = 0;
      str = 0;
      ch = -1;
      break;
    case 7: // pi0
      bary = 0;
      str = 0;
      ch = 0;
      break;
    case 11: // k+
      bary = 0;
      str = 1;
      ch = 1;
      break;
    case 13: // k-
      bary = 0;
      str = -1;
      ch = -1;
      break;
    case 15: // k0
      bary = 0;
      str = 1;
      ch = 0;
      break;
    case 17: // k0bar
      bary = 0;
      str = -1;
      ch = 0;
      break;
    case 21: // lambda
      bary = 1;
      str = -1;
      ch = 0;
      break;
    case 23: // sigma+
      bary = 1;
      str = -1;
      ch = 1;
      break;
    case 25: // sigma0
      bary = 1;
      str = -1;
      ch = 0;
      break;
    case 27: // sigma-
      bary = 1;
      str = -1;
      ch = -1;
      break;
    case 29: // xi0
      bary = 1;
      str = -2;
      ch = 0;
      break;
    case 31: // xi-
      bary = 1;
      str = -2;
      ch = -1;
      break;
    default:
      G4cout << " Unknown particle type " << type << G4endl;
    };

    Qnums[0] = bary;
    Qnums[1] = str;
    Qnums[2] = ch;
    return Qnums;
}



const G4double G4CascadeChannel::energyScale[31] = 
  { 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 1.0,  1.5, 2.0,  2.5, 3.0,  3.5, 4.0,  4.5, 5.0,
    5.5, 6.0,  6.5, 7.0,  7.5, 8.0,  8.5, 9.0,  9.5, 10.0, 15.0 };
