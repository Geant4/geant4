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
// $Id: G4PreCompoundFragment.hh 100691 2016-10-31 11:26:25Z gcosmo $
//
//  J. M. Quesada (August 2008).  
//  Based  on previous work by V. Lara
//
// Modified:
// 03.09.2008 by J. M. Quesada for external choice of inverse 
//                    cross section option (default OPTxs=2)
// 06.09.2008 by JMQ Also external choice has been added for
//     superimposed Coulomb barrier (if useSICB=true, default false) 
// 20.08.2010 V.Ivanchenko added int Z and A and cleanup; added 
//                        G4ParticleDefinition to constructor

#ifndef G4PreCompoundFragment_h
#define G4PreCompoundFragment_h 1

#include "G4VPreCompoundFragment.hh"

class G4PreCompoundFragment : public G4VPreCompoundFragment
{
public:  

  G4PreCompoundFragment(const G4ParticleDefinition*,
			G4VCoulombBarrier * aCoulombBarrier);
  
  virtual ~G4PreCompoundFragment();
      
  // ================================================
  // Methods for calculating the emission probability
  // ================================================
  
  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  G4double CalcEmissionProbability(const G4Fragment & aFragment);
  
  G4double SampleKineticEnergy(const G4Fragment & aFragment);

protected:

  virtual G4double GetAlpha() const = 0;

  virtual G4double GetBeta() const = 0;

  G4double CrossSection(G4double ekin) const;

  virtual G4double 
  ProbabilityDistributionFunction(G4double K, 
				  const G4Fragment & aFragment) = 0; 

private:	
  // This method performs integration for probability function over 
  // fragment kinetic energy
  G4double IntegrateEmissionProbability(G4double Low, G4double Up, 
					const G4Fragment & aFragment);	

  G4double GetOpt0(G4double ekin) const;

  // operators
  G4PreCompoundFragment(const G4PreCompoundFragment &right) = delete;
  const G4PreCompoundFragment& 
  operator= (const G4PreCompoundFragment &right) = delete;
  G4int operator==(const G4PreCompoundFragment &right) const = delete;
  G4int operator!=(const G4PreCompoundFragment &right) const = delete;

  G4int index;

  G4double muu;
  G4double probmax;
};

#endif
