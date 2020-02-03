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
//
// 
// Based on G4AtomicTransition.cc by 
// Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
// 
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// 4 Mar 2002: first implementation
//
// -------------------------------------------------------------------

#include "G4RDAugerTransition.hh"

// the final shell in wich the electron goes is needed, to know the data for the auger electron emitted
// (i.e. originating shell id, electron energy and transition probability) 

G4RDAugerTransition::G4RDAugerTransition(G4int finalShell, std::vector<G4int> transIds,
				     const std::map<G4int,std::vector<G4int>,std::less<G4int> >* idMap,
				     const std::map<G4int,G4DataVector,std::less<G4int> >* energyMap,
				     const std::map<G4int,G4DataVector,std::less<G4int> >* probabilityMap)
{
  finalShellId = finalShell;
  augerOriginatingShellIdsMap = *idMap;
  augerTransitionEnergiesMap = *energyMap;
  augerTransitionProbabilitiesMap = *probabilityMap;
  transitionOriginatingShellIds = transIds;
  

}

G4RDAugerTransition::~G4RDAugerTransition()
{ 

}

// Returns the ids of the shells from wich an auger electron culd came from, given th shell
// from wich the transition electron cames from.

const std::vector<G4int>* G4RDAugerTransition::AugerOriginatingShellIds(G4int startShellId) const
{
  std::map<G4int,std::vector<G4int>,std::less<G4int> >::const_iterator shellId = augerOriginatingShellIdsMap.find(startShellId);

  const std::vector<G4int>* dataSet = &(*shellId).second;
  //const std::vector<G4int>* dataOut = 0;

  if (dataSet->size() == 0) {G4cout << "Error: no auger Id found"<< G4endl;}
  else {
    
    // dataOut = &dataSet;  

  }

  return dataSet;
}

// Returns the ids of the shells from wich an electron cuuld fill the vacancy in finalShellId

const std::vector<G4int>* G4RDAugerTransition::TransitionOriginatingShellIds() const
{

  const std::vector<G4int>* dataSet = &transitionOriginatingShellIds;
  return dataSet;
}

// Returns the energiess of the possible auger electrons, given th shell
// from wich the transition electron cames from.

const G4DataVector* G4RDAugerTransition::AugerTransitionEnergies(G4int startShellId) const
{
  std::map<G4int,G4DataVector,std::less<G4int> >::const_iterator shellId = augerTransitionEnergiesMap.find(startShellId);
  const G4DataVector* dataSet = &(*shellId).second;


  return dataSet;
}

// Returns the emission probabilities of the auger electrons, given th shell
// from wich the transition electron cames from.

const G4DataVector* G4RDAugerTransition::AugerTransitionProbabilities(G4int startShellId) const
{
  std::map<G4int,G4DataVector,std::less<G4int> >::const_iterator shellId = augerTransitionProbabilitiesMap.find(startShellId);
  const G4DataVector* dataSet = &(*shellId).second;
  return dataSet; 
}

G4int G4RDAugerTransition::FinalShellId() const
{ 
  return finalShellId;
}

// Returns the id of the shell from wich come the auger electron , given the shell
// from wich the transition electron cames from and the index number.

G4int G4RDAugerTransition::AugerOriginatingShellId(G4int index, G4int startShellId) const
{
  const std::vector<G4int>* ids = AugerOriginatingShellIds(startShellId);
  // G4int i = 
  std::vector<G4int>::const_iterator pos = ids->begin();
  G4int n = 0;
  n = *(pos+index);
  return n;
}

// Returns the energy of the auger electron, given the shell
// from wich the transition electron cames from and the index number.

G4double G4RDAugerTransition::AugerTransitionEnergy(G4int index, G4int startShellId) const
{
  const G4DataVector* energies = AugerTransitionEnergies(startShellId);
    G4double energy = 0;    
  if (index < (G4int) energies->size()) {
    G4DataVector::const_iterator pos = energies->begin();
    energy = *(pos+index);
  }
  return energy; 
}

// Returns the probability of the auger emission, given the shell
// from wich the transition electron cames from and the index number.

G4double G4RDAugerTransition::AugerTransitionProbability(G4int index, G4int startShellId) const
{

  const G4DataVector *probabilities = AugerTransitionProbabilities(startShellId);
  G4DataVector::const_iterator pos = probabilities->begin();

  G4double probability = 0; 
  probability = *(pos+index);

  return  probability;
  
}

G4int G4RDAugerTransition::TransitionOriginatingShellId(G4int index) const
{
  return  transitionOriginatingShellIds[index];
}





















