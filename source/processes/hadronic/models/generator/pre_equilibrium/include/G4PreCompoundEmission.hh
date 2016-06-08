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
// $Id: G4PreCompoundEmission.hh,v 1.8 2001/12/13 12:04:18 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#ifndef G4PreCompoundEmission_h
#define G4PreCompoundEmission_h 1

#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundFragmentVector.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"
#include "G4PreCompoundParameters.hh"

class G4PreCompoundEmission
{
public:
  G4PreCompoundEmission(const G4Fragment& aFragment);
  ~G4PreCompoundEmission() {};

private:
  G4PreCompoundEmission() {};
  G4PreCompoundEmission(const G4PreCompoundEmission &right);
  const G4PreCompoundEmission& operator=(const G4PreCompoundEmission &right);
  G4bool operator==(const G4PreCompoundEmission &right) const;
  G4bool operator!=(const G4PreCompoundEmission &right) const;

public:

  void Initialize(const G4Fragment & aFragment) 
  {
    theFragmentsVector.Initialize(aFragment);
    return;
  }
	
  G4double GetTotalProbability(const G4Fragment & aFragment) 
  {
    return theFragmentsVector.CalculateProbabilities(aFragment);
  }
	
  G4ReactionProduct * PerformEmission(G4Fragment & aFragment);

	
private:

  G4ThreeVector AngularDistribution(G4VPreCompoundFragment * theFragment,
				    const G4Fragment& aFragment,
				    const G4double KineticEnergy) const;
		


  G4double rho(const G4double p, const G4double h, const G4double g, 
	       const G4double E, const G4double Ef) const;

  //  G4double bessi0(const G4double x) const;					 
					 
  // A vector with the allowed emission fragments 
  G4PreCompoundFragmentVector theFragmentsVector;
	
  // Projectile energy
  G4double ProjEnergy;

  // Projectile direction
  G4ThreeVector theIncidentDirection;


};
#endif
