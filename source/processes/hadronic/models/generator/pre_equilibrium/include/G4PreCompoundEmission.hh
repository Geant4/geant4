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
// $Id: G4PreCompoundEmission.hh,v 1.10 2003/03/24 13:56:44 larazb Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#ifndef G4PreCompoundEmission_h
#define G4PreCompoundEmission_h 1

#include "G4VPreCompoundFragment.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"
#include "G4PreCompoundParameters.hh"
#include "G4PreCompoundFragmentVector.hh"

class G4VPreCompoundEmissionFactory;


class G4PreCompoundEmission
{
public:
  G4PreCompoundEmission();
  ~G4PreCompoundEmission();

private:
  G4PreCompoundEmission(const G4PreCompoundEmission &right);
  const G4PreCompoundEmission& operator=(const G4PreCompoundEmission &right);
  G4bool operator==(const G4PreCompoundEmission &right) const;
  G4bool operator!=(const G4PreCompoundEmission &right) const;

public:

  inline void SetUp(const G4Fragment & aFragment);
  inline void Initialize(const G4Fragment & aFragment);
	
  inline G4double GetTotalProbability(const G4Fragment & aFragment);
	
  G4ReactionProduct * PerformEmission(G4Fragment & aFragment);

  void SetDefaultModel();
  void SetHETCModel();

private:

  G4ThreeVector AngularDistribution(G4VPreCompoundFragment * theFragment,
				    const G4Fragment& aFragment,
				    const G4double KineticEnergy) const;
		


  G4double rho(const G4double p, const G4double h, const G4double g, 
	       const G4double E, const G4double Ef) const;

  //==============
  // Data Members
  //==============

  // A vector with the allowed emission fragments 
  G4PreCompoundFragmentVector * theFragmentsVector;
  G4VPreCompoundEmissionFactory * theFragmentsFactory;

  // Projectile energy
  G4double ProjEnergy;

  // Projectile direction
  G4ThreeVector theIncidentDirection;

};

#include "G4PreCompoundEmission.icc"

#endif
