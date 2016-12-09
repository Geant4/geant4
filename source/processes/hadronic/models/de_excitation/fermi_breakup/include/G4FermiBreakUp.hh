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
// $Id: G4FermiBreakUp.hh 100379 2016-10-19 15:05:35Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4FermiBreakUp_h
#define G4FermiBreakUp_h 1

#include "G4VFermiBreakUp.hh"
#include "globals.hh"
#include "G4FermiConfiguration.hh"
#include "G4VFermiFragment.hh"
#include "G4FermiPhaseSpaceDecay.hh"
#include <vector>

class G4FermiFragmentsPool;
class G4Pow;

class G4FermiBreakUp : public G4VFermiBreakUp 
{
public:

  explicit G4FermiBreakUp();
  virtual ~G4FermiBreakUp();

  virtual void Initialise() final;

  virtual G4bool IsApplicable(G4int Z, G4int A, G4double mass) const final;

  // new interface - vector of products is added to the provided vector
  // primary fragment is deleted or is modified and added to the list
  // of products 
  virtual void BreakFragment(G4FragmentVector*, G4Fragment* theNucleus) final;
  
private:

  G4double CoulombBarrier(const std::vector<const G4VFermiFragment*>* v);

  G4double DecayProbability(G4int A, G4double TotalE, const G4FermiConfiguration*);

  const std::vector<const G4VFermiFragment*>*
  SelectConfiguration(G4int Z, G4int A, G4double mass);

  G4FermiBreakUp(const G4FermiBreakUp &right) = delete;  
  const G4FermiBreakUp & operator=(const G4FermiBreakUp &right) = delete;
  G4bool operator==(const G4FermiBreakUp &right) const = delete;
  G4bool operator!=(const G4FermiBreakUp &right) const = delete;

  G4FermiFragmentsPool* thePool;

  std::vector<G4double> NormalizedWeights;
  
  G4double Coef;
  G4double ConstCoeff;
  size_t   nmax;
  
  G4Pow* g4calc;

  const G4FermiPhaseSpaceDecay*        thePhaseSpace;
  std::vector<G4double>                massRes;
  std::vector<const G4VFermiFragment*> frag;  
};


#endif


