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
// $Id: G4PreCompoundAlpha.hh,v 1.9 2002/01/15 12:43:07 vlara Exp $
// GEANT4 tag $Name: geant4-04-00-patch-02 $
//
// by V. Lara 

#ifndef G4PreCompoundAlpha_h
#define G4PreCompoundAlpha_h 1

#include "G4VPreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Alpha.hh"

#include "G4AlphaCoulombBarrier.hh"

class G4PreCompoundAlpha : public G4VPreCompoundIon
{
public:
    // default constructor
    G4PreCompoundAlpha():G4VPreCompoundIon(4,2,&theAlphaCoulombBarrier,"alpha") {}

    // copy constructor
    G4PreCompoundAlpha(const G4PreCompoundAlpha &right): G4VPreCompoundIon(right) {}

    // destructor
    ~G4PreCompoundAlpha() {}

    // operators  
    const G4PreCompoundAlpha & operator=(const G4PreCompoundAlpha &right) {
	if (&right != this) this->G4VPreCompoundIon::operator=(right);
	return *this;
    }

    G4bool operator==(const G4PreCompoundAlpha &right) const
	{ return G4VPreCompoundIon::operator==(right);}
  
    G4bool operator!=(const G4PreCompoundAlpha &right) const
	{ return G4VPreCompoundIon::operator!=(right);}

    G4ReactionProduct * GetReactionProduct() const 
	{
	    G4ReactionProduct * theReactionProduct = 
		new G4ReactionProduct(G4Alpha::AlphaDefinition());
	    theReactionProduct->SetMomentum(GetMomentum().vect());
	    theReactionProduct->SetTotalEnergy(GetMomentum().e());
#ifdef pctest
	    theReactionProduct->SetCreatorModel("G4PrecompoundModel");
#endif
	    return theReactionProduct;
	}
    
public:
    void CalcExcitonLevelDensityRatios(const G4double Excitons,
				       const G4double Particles)
	{
	    // Level density ratios are calculated according to the formula
	    // (P!*(N-1)!)/((P-Af)!*(N-1-Af)!*Af!)
	    // where  P is number of particles
	    //        N is number of excitons
	    //        Af atomic number of emitting fragment
	    // the next is a simplification for alphas (Af = 4)
	    
	    SetExcitonLevelDensityRatio(((Particles*(Excitons-1.0))*
					 ((Particles-1.0)*(Excitons-2.0)/2.0)*
					 ((Particles-2.0)*(Excitons-3.0)/3.0)*
					 ((Particles-3.0)*(Excitons-4.0)/4.0))/6.0);
	}
    
    
    
    void CalcCondensationProbability(const G4double A)
	// This method computes condensation probability to create a fragment
	// consisting from N nucleons inside a nucleus with A nucleons 
	// This value comes from the formula N^3 (N/A)^(N-1) with N = 4 (alpha)
	{
	    SetCondensationProbability(4096.0/(A*A*A));
	}
    
private:
    
    virtual G4double GetCCoef(const G4double aZ) const;
    
    G4AlphaCoulombBarrier theAlphaCoulombBarrier;
    
};

#endif

inline G4double G4PreCompoundAlpha::GetCCoef(const G4double aZ) const
{
  G4double C = 0.0;

  if (aZ <= 30) {
    C = 0.10;
  } else if (aZ <= 50) {
    C = 0.1 + -((aZ-50.)/20.)*0.02;
  } else if (aZ < 70) {
    C = 0.08 + -((aZ-70.)/20.)*0.02;
  } else {
    C = 0.06;
  }
  return C;            
}

 
