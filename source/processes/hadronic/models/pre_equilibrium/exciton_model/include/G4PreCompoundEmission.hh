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
// $Id: G4PreCompoundEmission.hh,v 1.3 2006/06/29 20:58:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
