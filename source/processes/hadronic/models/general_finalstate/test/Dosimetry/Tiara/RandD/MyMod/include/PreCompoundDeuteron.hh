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
// $Id: PreCompoundDeuteron.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#ifndef PreCompoundDeuteron_h
#define PreCompoundDeuteron_h 1

#include "VPreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Deuteron.hh"

#include "G4DeuteronCoulombBarrier.hh"
#include "globals.hh"

class PreCompoundDeuteron : public VPreCompoundIon
{
public:
  // default constructor
  PreCompoundDeuteron():VPreCompoundIon(2,1,&theDeuteronCoulombBarrier,"Deuteron") {}

  // copy constructor
  PreCompoundDeuteron(const PreCompoundDeuteron &right): VPreCompoundIon(right) {}

  // destructor
  ~PreCompoundDeuteron() {}

  // operators  
  const PreCompoundDeuteron & operator=(const PreCompoundDeuteron &right) {
    if (&right != this) this->VPreCompoundIon::operator=(right);
    return *this;
  }

  G4bool operator==(const PreCompoundDeuteron &right) const
  { return VPreCompoundIon::operator==(right);}

  
  G4bool operator!=(const PreCompoundDeuteron &right) const
  { return VPreCompoundIon::operator!=(right);}


    G4ReactionProduct * GetReactionProduct() const
	{
            G4ReactionProduct * theReactionProduct =
                new G4ReactionProduct(G4Deuteron::DeuteronDefinition());
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
	    return 1.0 + C/2.0;
	}
    
    virtual G4double GetBeta(const G4Fragment& fragment)
	{
	    G4double Common = 1.5*fermi*pow(fragment.GetA()-2,1./3) +
1.89*fermi;
	    G4double z = fragment.GetZ()-1;
	    Common = -z*eplus*eplus*MeV*MeV/coulomb/coulomb/Common;
	    if(z<=10)
	      return Common*0.48*z*0.01;
	    else if(z<20)
	      return Common*(0.48 + (0.64-0.48)*(z-10)*0.1);
	    else if(z<30)
	      return Common*(0.64 + (0.74-0.64)*(z-20)*0.1);
	    else if(z<50)
	      return Common*(0.74 + (0.83-0.74)*(z-30)*0.05);
	    return Common*0.86;
	}

    virtual G4double FactorialFactor(const G4double N, const G4double P)
	{
	    return 
		(N-1.0)*(N-2.0)*(P-1.0)*P/2.0;
	}

    virtual G4double CoalescenceFactor(const G4double A)
	{
	    return 16.0/A;
	}    
private:

    G4DeuteronCoulombBarrier theDeuteronCoulombBarrier;

};

#endif
 

