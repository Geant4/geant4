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
// $Id: PreCompoundProton.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#ifndef PreCompoundProton_h
#define PreCompoundProton_h 1

#include "VPreCompoundNucleon.hh"
#include "G4ReactionProduct.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

#include "G4ProtonCoulombBarrier.hh"
#include "globals.hh"

class PreCompoundProton : public VPreCompoundNucleon
{
public:
  // default constructor
  PreCompoundProton():VPreCompoundNucleon(1,1,&theProtonCoulombBarrier,"Proton") {}

  // copy constructor
  PreCompoundProton(const PreCompoundProton &right): VPreCompoundNucleon(right) {}

  // destructor
  ~PreCompoundProton() {};

  // operators  
  const PreCompoundProton & operator=(const PreCompoundProton &right) {
    if (&right != this) this->VPreCompoundNucleon::operator=(right);
    return *this;
  };

  G4bool operator==(const PreCompoundProton &right) const
  { return VPreCompoundNucleon::operator==(right);}

  
  G4bool operator!=(const PreCompoundProton &right) const
  { return VPreCompoundNucleon::operator!=(right);}


    G4ReactionProduct * GetReactionProduct() const
	{
	    G4ReactionProduct * theReactionProduct = 
		new G4ReactionProduct(G4Proton::ProtonDefinition());
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
	    G4double aZ = G4double(GetRestZ());
	    G4double C = 0.0;
	    if (aZ >= 70) {
		C = 0.10;
	    } else {
		C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
	    }
	    return 1.0 + C;
	}

    virtual G4double GetBeta(const G4Fragment& fragment)
	{
	    G4double Common =
1.5*fermi*pow(fragment.GetA()-1,1./3)+0.88*fermi;
	    G4double z = fragment.GetZ()-1;
	    Common = -z*eplus*eplus*MeV*MeV/coulomb/coulomb/Common;
	    if(z<=10)
	      return Common*0.42*z*0.01;
	    else if(z<=20)
	      return Common*(0.42 + (0.58-0.42)*(z-10.)*0.1);
	    else if(z<=30)
	      return Common*(0.58 + (0.68-0.58)*(z-20.)*0.1);
	    else if(z<=50)
	      return Common*(0.68 + (0.77-0.68)*(z-30.)*0.05);
	    return Common*0.8;
	}

    virtual G4bool IsItPossible(const G4Fragment& aFragment)
	{
	    return (aFragment.GetNumberOfCharged() >= 1);
	}


private:

  G4ProtonCoulombBarrier theProtonCoulombBarrier;

};

#endif
 
