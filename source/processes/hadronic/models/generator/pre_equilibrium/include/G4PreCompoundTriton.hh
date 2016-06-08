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
// $Id: G4PreCompoundTriton.hh,v 1.9 2002/01/15 12:47:54 vlara Exp $
// GEANT4 tag $Name: geant4-04-00-patch-02 $
//
// by V. Lara

#ifndef G4PreCompoundTriton_h
#define G4PreCompoundTriton_h 1

#include "G4VPreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Triton.hh"

#include "G4TritonCoulombBarrier.hh"

#include "G4ProtonCoulombBarrier.hh"


class G4PreCompoundTriton : public G4VPreCompoundIon
{
public:
  // default constructor
  G4PreCompoundTriton():G4VPreCompoundIon(3,1,&theTritonCoulombBarrier,"Triton") {}

  // copy constructor
  G4PreCompoundTriton(const G4PreCompoundTriton &right): G4VPreCompoundIon(right) {}
	
  // destructor
  ~G4PreCompoundTriton() {}

  // operators  
  const G4PreCompoundTriton & operator=(const G4PreCompoundTriton &right) {
    if (&right != this) this->G4VPreCompoundIon::operator=(right);
    return *this;
  }

  G4bool operator==(const G4PreCompoundTriton &right) const
  { return G4VPreCompoundIon::operator==(right);}

  
  G4bool operator!=(const G4PreCompoundTriton &right) const
  { return G4VPreCompoundIon::operator!=(right);}

    G4ReactionProduct *  GetReactionProduct() const 
	{
	    G4ReactionProduct * theReactionProduct = 
		new G4ReactionProduct(G4Triton::TritonDefinition());
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
    // (P!*(N-1)!)/((P-Af)!*(N-1-Af)!*Af! (Af-1)!)
    // where  P is number of particles
    //        N is number of excitons
    //        Af atomic number of emitting fragment
    // the next is a simplification for tritons (Af = 3)

    SetExcitonLevelDensityRatio(((Particles*(Excitons-1.0))*
				 ((Particles-1.0)*(Excitons-2.0)/2.0)*
				 ((Particles-2.0)*(Excitons-3.0)/3.0)/2.0));
  }


  void CalcCondensationProbability(const G4double A)
    // This method computes condensation probability to create a fragment
    // consisting from N nucleons inside a nucleus with A nucleons 
    // This value comes from the formula N^3 (N/A)^(N-1) with N = 3 (triton)
  {
    SetCondensationProbability(243.0/(A*A));
  }


private:

  virtual G4double GetCCoef(const G4double aZ) const;

  G4TritonCoulombBarrier theTritonCoulombBarrier;
};


#endif


inline G4double G4PreCompoundTriton::GetCCoef(const G4double aZ) const
{
  G4double C = 0.0;

  if (aZ >= 70) {
    C = 0.10;
  } else {
    C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
  }

  return C/3.0;                               
}


