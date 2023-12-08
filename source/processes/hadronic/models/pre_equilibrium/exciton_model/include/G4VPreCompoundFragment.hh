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
//
// J. M. Quesada (August 2008).  
// Based  on previous work by V. Lara
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option 
// JMQ (06 September 2008) Also external choice has been added for:
//                      - superimposed Coulomb barrier (if useSICB=true) 
// 20.08.2010 V.Ivanchenko added int Z and A and cleanup; added 
//                        G4ParticleDefinition to constructor, 
//                        inline method to build G4ReactionProduct; 
//                        remove string name
//                        

#ifndef G4VPreCompoundFragment_h
#define G4VPreCompoundFragment_h 1

#include "G4ios.hh"
#include <iomanip>
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Fragment.hh"
#include "G4ReactionProduct.hh"
#include "G4Pow.hh"

class G4NuclearLevelData;
class G4DeexPrecoParameters;
class G4VCoulombBarrier;

class G4VPreCompoundFragment
{
public:  

  explicit G4VPreCompoundFragment(const G4ParticleDefinition*,
				  G4VCoulombBarrier*);
  
  virtual ~G4VPreCompoundFragment();
  
  friend std::ostream& 
  operator<<(std::ostream&, const G4VPreCompoundFragment*);
  friend std::ostream& 
  operator<<(std::ostream&, const G4VPreCompoundFragment&);
  
  // =====================
  // Pure Virtual methods
  // =====================
  
  // Initialization method
  void Initialize(const G4Fragment& aFragment);
    
  // Methods for calculating the emission probability
  // ------------------------------------------------
  
  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  virtual G4double CalcEmissionProbability(const G4Fragment&) = 0;
  
  // sample kinetic energy of emitted fragment
  virtual G4double SampleKineticEnergy(const G4Fragment&) = 0;

  inline G4bool IsItPossible(const G4Fragment& aFragment) const;
  
  inline G4ReactionProduct* GetReactionProduct() const; 	
  
  G4int GetA() const { return theA; }
  
  G4int GetZ() const { return theZ; }
  
  G4int GetRestA() const { return theResA; }
  
  G4int GetRestZ() const { return theResZ; }

  G4double GetBindingEnergy() const { return theBindingEnergy; }
  
  G4double GetEnergyThreshold() const
  {
    return theMaxKinEnergy - theCoulombBarrier;
  }

  G4double GetEmissionProbability() const { return theEmissionProbability; }
    
  G4double GetNuclearMass() const { return theMass; }
  
  G4double GetRestNuclearMass() const { return theResMass; }

  const G4LorentzVector& GetMomentum() const { return theMomentum; }
  
  void SetMomentum(const G4LorentzVector& lv) { theMomentum = lv; }
  
  //for inverse cross section choice
  void SetOPTxs(G4int opt) { OPTxs = opt; }
  //for superimposed Coulomb Barrier for inverse cross sections
  void UseSICB(G4bool use) { useSICB = use; } 

  G4VPreCompoundFragment(const G4VPreCompoundFragment &right) = delete;
  const G4VPreCompoundFragment& 
  operator= (const G4VPreCompoundFragment &right) = delete;  
  G4bool operator==(const G4VPreCompoundFragment &right) const = delete;
  G4bool operator!=(const G4VPreCompoundFragment &right) const = delete;

protected:

  virtual G4double GetAlpha() const = 0;

  virtual G4double GetBeta() const { return -theCoulombBarrier; }

  G4NuclearLevelData* fNucData;
  G4DeexPrecoParameters* theParameters;
  G4Pow* g4calc;

  G4int theA;
  G4int theZ;
  G4int theResA{0};
  G4int theResZ{0};
  G4int theFragA{0};
  G4int theFragZ{0};

  G4double theResA13{0.0};
  G4double theBindingEnergy{0.0};
  G4double theMinKinEnergy{0.0};
  G4double theMaxKinEnergy{0.0};
  G4double theResMass{0.0};
  G4double theReducedMass{0.0};
  G4double theMass;

  G4double theEmissionProbability{0.0};
  G4double theCoulombBarrier{0.0};

  //for inverse cross section choice
  G4int OPTxs{3};
  //for superimposed Coulomb Barrier for inverse cross sections
  G4bool useSICB{true};

private:

  const G4ParticleDefinition* particle;
  G4VCoulombBarrier* theCoulombBarrierPtr;
  G4LorentzVector theMomentum{0., 0., 0., 0.};
};

inline G4bool
G4VPreCompoundFragment::IsItPossible(const G4Fragment& aFragment) const
{
  G4int pplus = aFragment.GetNumberOfCharged();
  G4int pneut = aFragment.GetNumberOfParticles()-pplus;
  return (pneut >= theA - theZ && pplus >= theZ && theMaxKinEnergy > 0.0);
}

inline G4ReactionProduct* G4VPreCompoundFragment::GetReactionProduct() const
{
  G4ReactionProduct* theReactionProduct = new G4ReactionProduct(particle);
  theReactionProduct->SetMomentum(GetMomentum().vect());
  theReactionProduct->SetTotalEnergy(GetMomentum().e());
  return theReactionProduct;
}

#endif
