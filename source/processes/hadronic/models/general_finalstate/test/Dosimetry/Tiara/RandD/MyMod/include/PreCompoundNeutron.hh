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
//
// $Id: PreCompoundNeutron.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara


#ifndef PreCompoundNeutron_h
#define PreCompoundNeutron_h 1

#include "VPreCompoundNucleon.hh"
#include "G4ReactionProduct.hh"
#include "G4Neutron.hh"
#include "Randomize.hh"

#include "G4NeutronCoulombBarrier.hh"

class PreCompoundNeutron : public VPreCompoundNucleon
{
public:
  // default constructor
  PreCompoundNeutron() : VPreCompoundNucleon(1,0,&theNeutronCoulomBarrier,"Neutron") {}

  // copy constructor
  PreCompoundNeutron(const PreCompoundNeutron &right): VPreCompoundNucleon(right) {}

  // destructor
  ~PreCompoundNeutron() {}

  // operators  
  const PreCompoundNeutron & operator=(const PreCompoundNeutron &right) {
    if (&right != this) this->VPreCompoundNucleon::operator=(right);
    return *this;
  }

  G4bool operator==(const PreCompoundNeutron &right) const
  { return VPreCompoundNucleon::operator==(right);}
  
  G4bool operator!=(const PreCompoundNeutron &right) const
  { return VPreCompoundNucleon::operator!=(right);}


    G4ReactionProduct * GetReactionProduct() const
	{
	    G4ReactionProduct * theReactionProduct = 
		new G4ReactionProduct(G4Neutron::NeutronDefinition());
	    theReactionProduct->SetMomentum(GetMomentum().vect());
	    theReactionProduct->SetTotalEnergy(GetMomentum().e());
#ifdef pctest
	    theReactionProduct->SetCreatorModel("G4PrecompoundModel");
#endif
	    return theReactionProduct;
	}
    
private:
    virtual G4double GetAlpha()
	{
	    return 0.76+2.2/pow(GetRestA(),1.0/3.0);
	}

    virtual G4double GetBeta(const G4Fragment& fragment) 
	{
	    return (2.12/pow(GetRestA(),2.0/3.0)-0.05)*MeV/GetAlpha();
	}

    virtual G4bool IsItPossible(const G4Fragment& aFragment)
	{
	    return ((aFragment.GetNumberOfParticles()-aFragment.GetNumberOfCharged()) >= 1);  
	}
    

private:

  G4NeutronCoulombBarrier theNeutronCoulomBarrier;

};

#endif
 





