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
// $Id: PreCompoundHe3.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#ifndef PreCompoundHe3_h
#define PreCompoundHe3_h 1

#include "VPreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4He3.hh"

#include "G4He3CoulombBarrier.hh"
#include "globals.hh"

class PreCompoundHe3 : public VPreCompoundIon
{
public:
  // default constructor
  PreCompoundHe3():VPreCompoundIon(3,2,&theHe3CoulombBarrier,"He3") {}

  // copy constructor
  PreCompoundHe3(const PreCompoundHe3 &right): VPreCompoundIon(right) {}

  // destructor
  ~PreCompoundHe3() {}

  // operators  
  const PreCompoundHe3 & operator=(const PreCompoundHe3 &right) {
    if (&right != this) this->VPreCompoundIon::operator=(right);
    return *this;
  }

  G4bool operator==(const PreCompoundHe3 &right) const
  { return VPreCompoundIon::operator==(right);}

  
  G4bool operator!=(const PreCompoundHe3 &right) const
  { return VPreCompoundIon::operator!=(right);}


    G4ReactionProduct * GetReactionProduct() const
	{
            G4ReactionProduct * theReactionProduct =
                new G4ReactionProduct(G4He3::He3Definition());
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
	    if (aZ <= 30) {
		C = 0.10;
	    } else if (aZ <= 50) {
		C = 0.1 + -((aZ-50.)/20.)*0.02;
	    } else if (aZ < 70) {
		C = 0.08 + -((aZ-70.)/20.)*0.02;
	    } else {
		C = 0.06;
	    }
	    return 1.0 + C*(4.0/3.0);
	}

    virtual G4double GetBeta(const G4Fragment& fragment)
	{
	    G4double Common = 1.5*fermi*pow(fragment.GetA()-3,1./3.)+2.1*fermi;
	    G4double z = fragment.GetZ()-2;
	    Common = -2*z*eplus*eplus*MeV*MeV/coulomb/coulomb/Common;
	    if(z<=10)
	      return Common*0.62*z*0.1;
	    else if(z<=20)
	      return Common*(0.62+(0.76-0.62)*(z-10.)*0.1);
	    else if(z<=30)
	      return Common*(0.76 + (0.85-0.76)*(z-20.)*0.1);
	    else if(z<=50) return Common*(0.85 + (0.91-0.85)*(z-30.)*0.05);
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

    G4He3CoulombBarrier theHe3CoulombBarrier;

};

#endif
 

