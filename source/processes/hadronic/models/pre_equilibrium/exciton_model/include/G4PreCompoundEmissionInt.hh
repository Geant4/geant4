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
// GEANT4 Class header file
//
// File name:  G4PreCompoundEmissionInt
//
// Author:  V.Ivantchenko, 25 January 2025
//
// Class Description:
// Model implementation for pre-equilibrium emission of a light fragment
// from an excited nucleus. It is an alternative to the default model.
//

#ifndef G4PreCompoundEmissionInt_h
#define G4PreCompoundEmissionInt_h 1

#include "G4VPreCompoundFragment.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "G4PreCompoundFragmentVector.hh"

class G4VPreCompoundEmissionFactory;
class G4Pow;
class G4NuclearLevelData;

class G4PreCompoundEmissionInt
{
public:

  explicit G4PreCompoundEmissionInt(G4int verb);

  ~G4PreCompoundEmissionInt();

  void SetDefaultModel();

  void SetHETCModel();

  G4ReactionProduct* PerformEmission(G4Fragment& aFragment);

  inline G4double GetTotalProbability(const G4Fragment& aFragment);

  inline void SetOPTxs(G4int);

  inline void UseSICB(G4bool);

  G4PreCompoundEmissionInt(const G4PreCompoundEmissionInt& right) = delete;
  const G4PreCompoundEmissionInt& operator=
  (const G4PreCompoundEmissionInt& right) = delete;
  G4bool operator==(const G4PreCompoundEmissionInt& right) const = delete;
  G4bool operator!=(const G4PreCompoundEmissionInt& right) const = delete;

private:

  void AngularDistribution(G4VPreCompoundFragment* theFragment,
			   const G4Fragment& aFragment,
			   G4double kineticEnergy);
		
  G4double rho(G4int p, G4int h, G4double gg, 
               G4double E, G4double Ef) const;

  G4Pow* g4calc;
  G4NuclearLevelData* fNuclData;

  G4PreCompoundFragmentVector* theFragmentsVector{nullptr};
  G4VPreCompoundEmissionFactory* theFragmentsFactory{nullptr};

  // Momentum of emitted fragment
  G4ThreeVector theFinalMomentum;
  G4double fFermiEnergy;

  G4bool fUseAngularGenerator;

  G4int fModelID;
  G4int fVerbose;
};

inline G4double 
G4PreCompoundEmissionInt::GetTotalProbability(const G4Fragment& aFragment) 
{
  return theFragmentsVector->CalculateProbabilities(aFragment);
}

inline void G4PreCompoundEmissionInt::SetOPTxs(G4int opt)
{
  theFragmentsVector->SetOPTxs(opt);
}

inline void G4PreCompoundEmissionInt::UseSICB(G4bool use)
{
  theFragmentsVector->UseSICB(use);
}

#endif
