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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: PreCompoundEmission.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#ifndef PreCompoundEmission_h
#define PreCompoundEmission_h 1

#include "VPreCompoundFragment.hh"
#include "PreCompoundFragmentVector.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"
#include "PreCompoundParameters.hh"

class PreCompoundEmission
{
public:
  PreCompoundEmission();
  ~PreCompoundEmission() {};

private:
  PreCompoundEmission(const PreCompoundEmission &right);
  const PreCompoundEmission& operator=(const PreCompoundEmission &right);
  G4bool operator==(const PreCompoundEmission &right) const;
  G4bool operator!=(const PreCompoundEmission &right) const;

public:

  void Initialize(const G4Fragment & aFragment) 
  {
    ProjEnergy = aFragment.GetExcitationEnergy();
    theIncidentDirection = aFragment.GetMomentum().vect().unit();
    theFragmentsVector.Initialize(aFragment);
    return;
  }
	
  G4double GetTotalProbability(const G4Fragment & aFragment,G4double dLvlDensity) 
  {
    return theFragmentsVector.CalculateProbabilities(aFragment,dLvlDensity);
  }
	
  G4ReactionProduct * PerformEmission(G4Fragment & aFragment,G4double dLvlDensity);

	
private:

  G4ThreeVector AngularDistribution(VPreCompoundFragment * theFragment,
				    const G4Fragment& aFragment,
				    const G4double KineticEnergy,G4double dLvlDensity) const;
		


  G4double rho(const G4double p, const G4double h, const G4double g, 
	       const G4double E, const G4double Ef) const;

  //  G4double bessi0(const G4double x) const;					 
					 
  // A vector with the allowed emission fragments 
  PreCompoundFragmentVector theFragmentsVector;
	
  // Projectile energy
  G4double ProjEnergy;

  // Projectile direction
  G4ThreeVector theIncidentDirection;


};
#endif
