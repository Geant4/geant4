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
// $Id: PreCompoundTriton.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#ifndef PreCompoundTriton_h
#define PreCompoundTriton_h 1

#include "VPreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Triton.hh"

#include "G4TritonCoulombBarrier.hh"


class PreCompoundTriton : public VPreCompoundIon
{
public:
  // default constructor
  PreCompoundTriton():VPreCompoundIon(3,1,&theTritonCoulombBarrier,"Triton") {}

  // copy constructor
  PreCompoundTriton(const PreCompoundTriton &right): VPreCompoundIon(right) {}

  // destructor
  ~PreCompoundTriton() {}

  // operators  
  const PreCompoundTriton & operator=(const PreCompoundTriton &right) {
    if (&right != this) this->VPreCompoundIon::operator=(right);
    return *this;
  }

  G4bool operator==(const PreCompoundTriton &right) const
  { return VPreCompoundIon::operator==(right);}

  
  G4bool operator!=(const PreCompoundTriton &right) const
  { return VPreCompoundIon::operator!=(right);}


    G4ReactionProduct * GetReactionProduct() const
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
    
private:
    virtual G4double GetAlpha()
	{
	    G4double C = 0.0;
	    G4double aZ = GetZ() + GetRestZ();
	    if (aZ >= 70) {
		C = 0.10;
	    } else {
		C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375; 
	    }
 
	    return 1.0 + C/3.0;
	}

    virtual G4double GetBeta(const G4Fragment& fragment)
	{
	    G4double Common = 1.5*fermi*pow(fragment.GetA()-3,1./3)+2.1*fermi;
	    G4double z = fragment.GetZ()-1;
	    Common = -z*eplus*eplus*MeV*MeV/coulomb/coulomb/Common;
	    if(z<=10)
	      return Common*0.54*z*0.1;
	    else if(z<=20)
	      return Common*(0.54 + (0.70-0.54)*(z-10)*0.1);
	    else if(z<=30)
	      return Common*(0.70 + (0.80-0.70)*(z-20)*0.1);
	    else if(z<50)
	      return Common*(0.80 + (0.89-0.80)*(z-30)*0.05);
	    return Common*0.92;
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
 

