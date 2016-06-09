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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// 23.04.2011 V.Ivanchenko: make this class to be responsible for
//            selection of decay channel and decay
//

#ifndef G4FermiConfigurationList_h
#define G4FermiConfigurationList_h 1

#include "globals.hh"
#include "G4FermiConfiguration.hh"
#include "G4VFermiFragment.hh"
#include "G4FermiPhaseSpaceDecay.hh"
#include <vector>

class G4FermiFragmentsPool;
class G4Pow;

class G4FermiConfigurationList 
{
public:

  G4FermiConfigurationList();

  ~G4FermiConfigurationList();

  G4FragmentVector* GetFragments(const G4Fragment& theNucleus);
  
private:

  G4double CoulombBarrier(const std::vector<const G4VFermiFragment*>& v);

  G4double DecayProbability(G4int A, G4double TotalE, G4FermiConfiguration*);

  const std::vector<const G4VFermiFragment*>*
  SelectConfiguration(G4int Z, G4int A, G4double mass);

  G4FermiConfigurationList(const G4FermiConfigurationList &right);  
  const G4FermiConfigurationList & operator=(const G4FermiConfigurationList &right);
  G4bool operator==(const G4FermiConfigurationList &right) const;
  G4bool operator!=(const G4FermiConfigurationList &right) const;
  
  G4FermiFragmentsPool* thePool;

  std::vector<G4double> NormalizedWeights;
  
  // Kappa = V/V_0 it is used in calculation of Coulomb energy
  static const G4double Kappa;
  
  // Nuclear radius r0 (is a model parameter)
  static const G4double r0;

  G4double Coef;
  G4double ConstCoeff;
  size_t   nmax;
  
  G4Pow* g4pow;

  G4FermiPhaseSpaceDecay thePhaseSpace;

};


#endif


