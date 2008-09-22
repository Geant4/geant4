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
//J. M. Quesada (August 2008).  
//Based  on previous work by V. Lara
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option (default OPTxs=2)
// JMQ (06 September 2008) Also external choice has been added for
// superimposed Coulomb barrier (if useSICB=true, default false) 


#ifndef G4PreCompoundFragment_h
#define G4PreCompoundFragment_h 1

#include "G4VPreCompoundFragment.hh"

class G4PreCompoundFragment : public G4VPreCompoundFragment
{
protected:
  // default constructor
  G4PreCompoundFragment() {};
    
public:  
  // copy constructor
  G4PreCompoundFragment(const G4PreCompoundFragment &right);
    
  // constructor  
  G4PreCompoundFragment(const G4double anA, const G4double aZ,
			G4VCoulombBarrier * aCoulombBarrier,
			const G4String &  aName);
  
  virtual ~G4PreCompoundFragment();
  
  // ==========
  // operators 
  // ========== 
  
  const G4PreCompoundFragment& 
  operator= (const G4PreCompoundFragment &right);
  
  G4int operator==(const G4PreCompoundFragment &right) const;
  
  G4int operator!=(const G4PreCompoundFragment &right) const;
      
public:
  
  // Initialization method
//  void Initialize(const G4Fragment & aFragment);
    
  // ================================================
  // Methods for calculating the emission probability
  // ================================================
  
  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  G4double CalcEmissionProbability(const G4Fragment & aFragment);
  
  G4double GetKineticEnergy(const G4Fragment & aFragment);

private:	
  // This method performs integration for probability function over 
  // fragment kinetic energy
  G4double IntegrateEmissionProbability(const G4double & Low, 
					const G4double & Up, 
					const G4Fragment & aFragment);	
protected:

  virtual G4double 
  ProbabilityDistributionFunction(const G4double K, 
				  const G4Fragment & aFragment) = 0; 

};

#endif
