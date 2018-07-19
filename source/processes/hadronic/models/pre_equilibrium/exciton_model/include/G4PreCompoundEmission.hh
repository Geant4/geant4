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
// $Id: G4PreCompoundEmission.hh 107062 2017-11-01 15:01:02Z gcosmo $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choice has been added for:
//                      - superimposed Coulomb barrier (if useSICB=true) 
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup; 
//                         inline methods moved from icc file to hh

#ifndef G4PreCompoundEmission_h
#define G4PreCompoundEmission_h 1

#include "G4VPreCompoundFragment.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "G4PreCompoundFragmentVector.hh"

class G4VPreCompoundEmissionFactory;
class G4Pow;

class G4PreCompoundEmission
{
public:

  G4PreCompoundEmission();

  ~G4PreCompoundEmission();

  void SetDefaultModel();

  void SetHETCModel();

  G4ReactionProduct* PerformEmission(G4Fragment & aFragment);

  inline G4double GetTotalProbability(const G4Fragment & aFragment);

  inline void SetOPTxs(G4int);

  inline void UseSICB(G4bool);

private:

  void AngularDistribution(G4VPreCompoundFragment * theFragment,
			   const G4Fragment& aFragment,
			   G4double KineticEnergy);
		
  G4double rho(G4int p, G4int h, G4double gg, 
	       G4double E, G4double Ef) const;

  G4PreCompoundEmission(const G4PreCompoundEmission &right);
  const G4PreCompoundEmission& operator=(const G4PreCompoundEmission &right);
  G4bool operator==(const G4PreCompoundEmission &right) const;
  G4bool operator!=(const G4PreCompoundEmission &right) const;

  //==============
  // Data Members
  //==============

  G4Pow* g4calc;
  G4double fLevelDensity;
  G4double fFermiEnergy;

  // A vector with the allowed emission fragments 
  G4PreCompoundFragmentVector * theFragmentsVector;
  G4VPreCompoundEmissionFactory * theFragmentsFactory;

  // Momentum of emitted fragment
  G4ThreeVector theFinalMomentum;
  G4bool fUseAngularGenerator;
};

inline G4double 
G4PreCompoundEmission::GetTotalProbability(const G4Fragment& aFragment) 
{
  return theFragmentsVector->CalculateProbabilities(aFragment);
}

inline void G4PreCompoundEmission::SetOPTxs(G4int opt)
{
  theFragmentsVector->SetOPTxs(opt);
}

inline void G4PreCompoundEmission::UseSICB(G4bool use)
{
  theFragmentsVector->UseSICB(use);
}

#endif
