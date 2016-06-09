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
// by V. Lara

#ifndef G4HETCTriton_h
#define G4HETCTriton_h 1

#include "G4HETCChargedFragment.hh"
#include "G4ReactionProduct.hh"
#include "G4Triton.hh"

#include "G4TritonCoulombBarrier.hh"


class G4HETCTriton : public G4HETCChargedFragment
{
public:
  // default constructor
  G4HETCTriton():G4HETCChargedFragment(3,1,&theTritonCoulombBarrier,"Triton") {}

  // copy constructor
  G4HETCTriton(const G4HETCTriton &right): G4HETCChargedFragment(right) {}

  // destructor
  ~G4HETCTriton() {}

  // operators  
  const G4HETCTriton & operator=(const G4HETCTriton &right) 
  {
    if (&right != this) this->G4HETCChargedFragment::operator=(right);
    return *this;
  }

  G4bool operator==(const G4HETCTriton &right) const
  { 
    return G4HETCChargedFragment::operator==(right);
  }

  
  G4bool operator!=(const G4HETCTriton &right) const
  { 
    return G4HETCChargedFragment::operator!=(right);
  }


  G4ReactionProduct * GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct =
      new G4ReactionProduct(G4Triton::TritonDefinition());
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
    G4double C = 0.0;
    G4double aZ = GetZ() + GetRestZ();
    if (aZ >= 70) 
      {
	C = 0.10;
      } 
    else 
      {
	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375; 
      }
    
    return 1.0 + C/3.0;
  }
  
  virtual G4double GetBeta()
  {
    return -GetCoulombBarrier();
  }



  virtual G4double GetSpinFactor()
  {
    // 2s+1
    return 2.0;
  }

  virtual G4double K(const G4Fragment& aFragment);

private:
  
  G4TritonCoulombBarrier theTritonCoulombBarrier;
  
};

#endif
 

