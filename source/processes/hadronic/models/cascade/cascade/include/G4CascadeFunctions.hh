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
#ifndef G4_CASCADE_FUNCTIONS_HH
#define G4_CASCADE_FUNCTIONS_HH

#include <vector>
#include "globals.hh"
#include "G4CascadeChannel.hh"

template <class T>
class G4CascadeFunctions
{
public:
  static G4double getCrossSection(double ke);
  static G4int getMultiplicity(G4double ke);
  static std::vector<G4int> getOutgoingParticleTypes(G4int mult, G4double ke);
};

template <class T>
inline
G4double 
G4CascadeFunctions<T>::getCrossSection(double ke)
{
  std::pair<G4int, G4double> epair = G4CascadeChannel::interpolateEnergy(ke);
  G4int k = epair.first;
  G4double fraction = epair.second;
  return T::data.tot[k] + fraction*(T::data.tot[k+1] - T::data.tot[k]);
}

template <class T>
inline
G4int 
G4CascadeFunctions<T>::getMultiplicity(G4double ke)
{
  G4double multint(0.0);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = G4CascadeChannel::interpolateEnergy(ke);
  G4int k = epair.first;
  G4double fraction = epair.second;

  for (G4int m = 0; m < 6; ++m)
    {
      multint = T::data.multiplicities[m][k]
	+ fraction * (T::data.multiplicities[m][k+1] - T::data.multiplicities[m][k]);
      sigma.push_back(multint);
    }
  
  return G4CascadeChannel::sampleFlat(sigma);
}

template <class T>
inline
std::vector<G4int> 
G4CascadeFunctions<T>::getOutgoingParticleTypes(G4int mult, G4double ke)
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;
  
  std::pair<G4int, G4double> epair = G4CascadeChannel::interpolateEnergy(ke);
  G4int k = epair.first;
  G4double fraction = epair.second;

  G4int start = T::data.index[mult-2][0];
  G4int stop = T::data.index[mult-2][1];
 
  for (i = start; i < stop; i++) {
    sigint = T::data.crossSections[i][k] 
      + fraction*(T::data.crossSections[i][k+1] - T::data.crossSections[i][k]);
    sigma.push_back(sigint);
  }
 
  G4int channel = G4CascadeChannel::sampleFlat(sigma);

  std::vector<G4int> kinds;

  if (mult == 2) {
    for(i = 0; i < mult; i++) kinds.push_back(T::data.x2bfs[channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T::data.x3bfs[channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T::data.x4bfs[channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T::data.x5bfs[channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T::data.x6bfs[channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T::data.x7bfs[channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }

  return kinds;
}


#endif
