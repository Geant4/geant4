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
// $Id: G4PreCompoundFragmentVector.hh 100378 2016-10-19 15:03:27Z gcosmo $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara
//
// Modified:
// 03.09.2008 by J. M. Quesada for external choice of inverse 
// cross section option 
// 06.09.2008 JMQ Also external choice has been added for:
//                - superimposed Coulomb barrier (if useSICB=true) 
// 27.08.2010 V.Ivanchenko simplify and make more efficient by adding extra
//            vector of probabilities, moved constructor and destructor to source, 
//            simplify run time computations making inlined
// 

#ifndef G4PreCompoundFragmentVector_h
#define G4PreCompoundFragmentVector_h 1

#include "G4VPreCompoundFragment.hh"
#include "G4DataVector.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <vector>

typedef std::vector<G4VPreCompoundFragment*>  pcfvector;

class G4PreCompoundFragmentVector 
{
public:

  explicit G4PreCompoundFragmentVector(pcfvector * avector);

  ~G4PreCompoundFragmentVector();

  void SetVector(pcfvector * avector);

  void SetOPTxs(G4int);

  void UseSICB(G4bool);

  inline G4double CalculateProbabilities(const G4Fragment & aFragment);
	
  inline G4VPreCompoundFragment * ChooseFragment();
		  
private:

  G4PreCompoundFragmentVector(const G4PreCompoundFragmentVector &right) = delete;
  const G4PreCompoundFragmentVector& 
  operator=(const G4PreCompoundFragmentVector &right) = delete;
  G4bool operator==(const G4PreCompoundFragmentVector &right) const = delete;
  G4bool operator!=(const G4PreCompoundFragmentVector &right) const = delete;
  
  pcfvector * theChannels;
  G4DataVector probabilities;

  G4int nChannels;
};

inline G4double 
G4PreCompoundFragmentVector::CalculateProbabilities(const G4Fragment & aFragment)
{
  //G4cout << "## G4PreCompoundFragmentVector::CalculateProbabilities" << G4endl;
  G4double probtot = 0.0; 
  G4double prob; 
  for (G4int i=0; i< nChannels; ++i) { 
    (*theChannels)[i]->Initialize(aFragment);
    prob = 0.0;
    if ((*theChannels)[i]->IsItPossible(aFragment)) {
      prob = (*theChannels)[i]->CalcEmissionProbability(aFragment);
    }
    probtot += prob;
    probabilities[i] = probtot;
    //G4cout<<" prob= "<<prob<<" for "<<(*theChannels)[i]->GetName()<<G4endl;
  }
  return probtot;
}

inline G4VPreCompoundFragment* G4PreCompoundFragmentVector::ChooseFragment()
{
  G4double x = probabilities[nChannels-1]*G4UniformRand();
  G4int i=0;
  for (; i<nChannels; ++i) { 
    if(x <= probabilities[i]) { break; }
  }
  return (*theChannels)[i];  
}

#endif

