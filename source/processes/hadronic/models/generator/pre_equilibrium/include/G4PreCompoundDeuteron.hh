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
// $Id: G4PreCompoundDeuteron.hh,v 1.6.2.1 2001/06/28 19:13:31 gunter Exp $
// GEANT4 tag $Name:  $
//
// by V. Lara 

#ifndef G4PreCompoundDeuteron_h
#define G4PreCompoundDeuteron_h 1

#include "G4VPreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Deuteron.hh"

#include "G4DeuteronCoulombBarrier.hh"

class G4PreCompoundDeuteron : public G4VPreCompoundIon
{
public:
  // default constructor
  G4PreCompoundDeuteron():G4VPreCompoundIon(2,1,&theDeuteronCoulombBarrier,"Deuteron") {};
  // copy constructor
  G4PreCompoundDeuteron(const G4PreCompoundDeuteron &right):
    G4VPreCompoundIon(right) {}
	
  // destructor
  ~G4PreCompoundDeuteron() {}

  // operators  
  const G4PreCompoundDeuteron & operator=(const G4PreCompoundDeuteron &right) {
    if (&right != this) this->G4VPreCompoundIon::operator=(right);
    return *this;
  }

  G4bool operator==(const G4PreCompoundDeuteron &right) const
  { return G4VPreCompoundIon::operator==(right);}
  
  G4bool operator!=(const G4PreCompoundDeuteron &right) const
  { return G4VPreCompoundIon::operator!=(right);}

  G4ReactionProduct * GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct = new G4ReactionProduct(G4Deuteron::DeuteronDefinition());
    theReactionProduct->SetMomentum(GetMomentum().vect());
    theReactionProduct->SetTotalEnergy(GetMomentum().e());
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
    // the next is a simplification for deuterons (Af = 2)

    SetExcitonLevelDensityRatio(Particles*(Excitons-1.0)*
				(Particles-1.0)*(Excitons-2.0)/2.0);
  }


  void CalcCondensationProbability(const G4double A)
    // This method computes condensation probability to create a fragment
    // consisting from N nucleons inside a nucleus with A nucleons 
    // This value comes from the formula N^3 (N/A)^(N-1) with N = 2 (deuteron)
  {
    SetCondensationProbability(16.0/A);
  }

private:

  virtual G4double GetCCoef(const G4double aZ) const;

  G4DeuteronCoulombBarrier theDeuteronCoulombBarrier;

};

#endif

inline G4double G4PreCompoundDeuteron::GetCCoef(const G4double aZ) const
{
  G4double C = 0.0;

  if (aZ >= 70) {
    C = 0.10;
  } else {
    C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
  }
  return C/2.0;
}


 
