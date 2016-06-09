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
// $Id: G4PreCompoundTriton.hh,v 1.1 2003/08/26 18:54:21 lara Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// by V. Lara

#ifndef G4PreCompoundTriton_h
#define G4PreCompoundTriton_h 1

#include "G4PreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Triton.hh"

#include "G4TritonCoulombBarrier.hh"


class G4PreCompoundTriton : public G4PreCompoundIon
{
public:
  // default constructor
  G4PreCompoundTriton():G4PreCompoundIon(3,1,&theTritonCoulombBarrier,"Triton") {}

  // copy constructor
  G4PreCompoundTriton(const G4PreCompoundTriton &right): G4PreCompoundIon(right) {}

  // destructor
  ~G4PreCompoundTriton() {}

  // operators  
  const G4PreCompoundTriton & operator=(const G4PreCompoundTriton &right) {
    if (&right != this) this->G4PreCompoundIon::operator=(right);
    return *this;
  }

  G4bool operator==(const G4PreCompoundTriton &right) const
  { return G4PreCompoundIon::operator==(right);}

  
  G4bool operator!=(const G4PreCompoundTriton &right) const
  { return G4PreCompoundIon::operator!=(right);}


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

  virtual G4double FactorialFactor(const G4double N, const G4double P)
  {
    return 
      (N-3.0)*(P-2.0)*(
		       (((N-2.0)*(P-1.0))/2.0) *(
						 (((N-1.0)*P)/3.0) 
						 )
		       );
  }

  virtual G4double CoalescenceFactor(const G4double A)
  {
    return 243.0/(A*A);
  }    
private:

  G4TritonCoulombBarrier theTritonCoulombBarrier;

};

#endif
