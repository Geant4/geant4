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
//
// $Id: G4PreCompoundProton.hh,v 1.1 2003/08/26 18:54:18 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// by V. Lara

#ifndef G4PreCompoundProton_h
#define G4PreCompoundProton_h 1

#include "G4PreCompoundNucleon.hh"
#include "G4ReactionProduct.hh"
#include "G4Proton.hh"
#include "G4PreCompoundParameters.hh"
#include "Randomize.hh"

#include "G4ProtonCoulombBarrier.hh"


class G4PreCompoundProton : public G4PreCompoundNucleon
{
public:
  // default constructor
  G4PreCompoundProton():G4PreCompoundNucleon(1,1,&theProtonCoulombBarrier,"Proton") {}

  // copy constructor
  G4PreCompoundProton(const G4PreCompoundProton &right): G4PreCompoundNucleon(right) {}

  // destructor
  ~G4PreCompoundProton() {};

  // operators  
  const G4PreCompoundProton & operator=(const G4PreCompoundProton &right) {
    if (&right != this) this->G4PreCompoundNucleon::operator=(right);
    return *this;
  };

  G4bool operator==(const G4PreCompoundProton &right) const
  { return G4PreCompoundNucleon::operator==(right);}

  
  G4bool operator!=(const G4PreCompoundProton &right) const
  { return G4PreCompoundNucleon::operator!=(right);}


  G4ReactionProduct * GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct = 
      new G4ReactionProduct(G4Proton::ProtonDefinition());
    theReactionProduct->SetMomentum(GetMomentum().vect());
    theReactionProduct->SetTotalEnergy(GetMomentum().e());
#ifdef PRECOMPOUND_TEST
    theReactionProduct->SetCreatorModel("G4PrecompoundModel");
#endif
    return theReactionProduct;
  }
  
private:
  virtual G4double GetAlpha()
  {
    G4double aZ = static_cast<G4double>(GetRestZ());
    G4double C = 0.0;
    if (aZ >= 70) 
      {
	C = 0.10;
      } 
    else 
      {
	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
      }
    return 1.0 + C;
  }

  virtual G4double GetBeta()
  {
    return -GetCoulombBarrier();
  }
  
  virtual G4bool IsItPossible(const G4Fragment& aFragment)
  {
    return (aFragment.GetNumberOfCharged() >= 1);
  }
  
  
private:
  
  G4ProtonCoulombBarrier theProtonCoulombBarrier;
  
};

#endif
 
