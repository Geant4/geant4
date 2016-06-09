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


#ifndef G4HETCNeutron_h
#define G4HETCNeutron_h 1

#include "G4HETCFragment.hh"
#include "G4ReactionProduct.hh"
#include "G4Neutron.hh"
#include "G4PreCompoundParameters.hh"
#include "Randomize.hh"

#include "G4NeutronCoulombBarrier.hh"


class G4HETCNeutron : public G4HETCFragment
{
public:
  // default constructor
  G4HETCNeutron() : G4HETCFragment(1,0,&theNeutronCoulomBarrier,"Neutron") {}

  // copy constructor
  G4HETCNeutron(const G4HETCNeutron &right): G4HETCFragment(right) {}

  // destructor
  ~G4HETCNeutron() {}

  // operators  
  const G4HETCNeutron & operator=(const G4HETCNeutron &right) {
    if (&right != this) this->G4HETCFragment::operator=(right);
    return *this;
  }

  G4bool operator==(const G4HETCNeutron &right) const
  { return G4HETCFragment::operator==(right);}
  
  G4bool operator!=(const G4HETCNeutron &right) const
  { return G4HETCFragment::operator!=(right);}


  G4ReactionProduct * GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct = 
      new G4ReactionProduct(G4Neutron::NeutronDefinition());
    theReactionProduct->SetMomentum(GetMomentum().vect());
    theReactionProduct->SetTotalEnergy(GetMomentum().e());
#ifdef PRECOMPOUND_TEST
    theReactionProduct->SetCreatorModel("G4PrecompoundModel");
#endif
    return theReactionProduct;
  }

  virtual G4double GetKineticEnergy(const G4Fragment & aFragment);
  
private:
  virtual G4double GetAlpha()
  {
    return 0.76+2.2/pow(GetRestA(),1.0/3.0);
  }
  
  virtual G4double GetBeta() 
  {
    return (2.12/pow(GetRestA(),2.0/3.0)-0.05)*MeV/GetAlpha();
  }

  virtual G4double GetSpinFactor()
  {
    // (2s+1)
    return 2.0;
  }
  
  virtual G4double K(const G4Fragment& aFragment);

private:
  
  G4NeutronCoulombBarrier theNeutronCoulomBarrier;

};

#endif
 





