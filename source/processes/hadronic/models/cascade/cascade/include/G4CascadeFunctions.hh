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
// $Id: G4CascadeFunctions.hh,v 1.3 2010-05-14 18:28:02 mkelsey Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100407  M. Kelsey -- Return particle types std::vector<> by const ref,
//		using a static variable in the function as a buffer.
// 20100505  M. Kelsey -- Use new interpolator class, drop std::pair<>, move
//		sampleFlat(...) from G4CascadeChannel, move functionality
//		to new base class, to allow data-member buffers.  Move
//		function definitions to .icc file (needed with templating).

#ifndef G4_CASCADE_FUNCTIONS_HH
#define G4_CASCADE_FUNCTIONS_HH

#include <vector>
#include "globals.hh"
#include "G4CascadeSampler.hh"

template <class T>
class G4CascadeFunctions : public G4CascadeSampler {
public:
  static G4double getCrossSection(double ke) {
    return instance.findCrossSection(ke, T::data.tot);
  }

  static G4int getMultiplicity(G4double ke) {
    return instance.findMultiplicity(ke, T::data.multiplicities);
  }

  static const std::vector<G4int>& 
  getOutgoingParticleTypes(G4int mult, G4double ke);

private:
  G4CascadeFunctions() : G4CascadeSampler() {}
  static const G4CascadeFunctions<T> instance;
};

// Make sure singleton is instantiated
template <class T>
const G4CascadeFunctions<T> G4CascadeFunctions<T>::instance;


template <class T> inline 
const std::vector<G4int>& 
G4CascadeFunctions<T>::getOutgoingParticleTypes(G4int mult, G4double ke) {
  G4int channel = instance.findFinalStateIndex(mult, ke, T::data.index,
					       T::data.crossSections);

  static std::vector<G4int> kinds(8);	// FIXME:  This is not thread-safe!
  kinds.clear();

  G4int i;
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

#endif	/* G4_CASCADE_FUNCTIONS_HH */
