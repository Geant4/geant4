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

  ~G4PreCompoundFragmentVector() = default;

  void SetVector(pcfvector * avector);

  void SetOPTxs(G4int);

  void UseSICB(G4bool);

  G4double CalculateProbabilities(const G4Fragment & aFragment);
	
  G4VPreCompoundFragment * ChooseFragment();
		  
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

#endif

