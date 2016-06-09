//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4PreCompoundAlpha.hh,v 1.3 2006/06/29 20:58:06 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// by V. Lara

#ifndef G4PreCompoundAlpha_h
#define G4PreCompoundAlpha_h 1

#include "G4PreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Alpha.hh"

#include "G4AlphaCoulombBarrier.hh"


class G4PreCompoundAlpha : public G4PreCompoundIon
{
public:
  // default constructor
  G4PreCompoundAlpha():G4PreCompoundIon(4,2,&theAlphaCoulombBarrier,"Alpha") {}

  // copy constructor
  G4PreCompoundAlpha(const G4PreCompoundAlpha &right): G4PreCompoundIon(right) {}

  // destructor
  ~G4PreCompoundAlpha() {}

  // operators  
  const G4PreCompoundAlpha & operator=(const G4PreCompoundAlpha &right) {
    if (&right != this) this->G4PreCompoundIon::operator=(right);
    return *this;
  };

  G4bool operator==(const G4PreCompoundAlpha &right) const
  { return G4PreCompoundIon::operator==(right);}

  
  G4bool operator!=(const G4PreCompoundAlpha &right) const
  { return G4PreCompoundIon::operator!=(right);}


  G4ReactionProduct * GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct =
      new G4ReactionProduct(G4Alpha::AlphaDefinition());
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
    if (aZ <= 30) 
      {
	C = 0.10;
      } 
    else if (aZ <= 50) 
      {
	C = 0.1 + -((aZ-50.)/20.)*0.02;
      } 
    else if (aZ < 70) 
      {
	C = 0.08 + -((aZ-70.)/20.)*0.02;
      }
    else 
      {
	C = 0.06;
      }
    return 1.0+C;
  }

  virtual G4double GetBeta()
  {
    return -GetCoulombBarrier();
  }

  virtual G4double FactorialFactor(const G4double N, const G4double P)
  {
    return 
      (N-4.0)*(P-3.0)*(
		       (((N-3.0)*(P-2.0))/2.0) *(
						 (((N-2.0)*(P-1.0))/3.0) *(
									   (((N-1.0)*P)/2.0)
									   )
						 )
		       );
  }

  virtual G4double CoalescenceFactor(const G4double A)
  {
    return 4096.0/(A*A*A);
  }    
private:

  G4AlphaCoulombBarrier theAlphaCoulombBarrier;

};

#endif
 

