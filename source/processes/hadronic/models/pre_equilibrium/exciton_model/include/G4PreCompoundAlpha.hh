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
// $Id: G4PreCompoundAlpha.hh,v 1.9 2008-07-23 18:25:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara
//
//J. M. Quesada (Dic. 07) Added combinatorial factor Rj. Removed Coulomb barrier , GetAlpha and GetBeta methods


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

//JMQ (Sep. 07) combinatorial factor Rj
  virtual G4double GetRj(const G4int NumberParticles, const G4int NumberCharged)
  {
    G4double rj = 0.0;
    G4double denominator = NumberParticles*(NumberParticles-1)*(NumberParticles-2)*(NumberParticles-3);
    if(NumberCharged >=2 && (NumberParticles-NumberCharged) >=2 ) rj = 6.0*static_cast<G4double>(NumberCharged*(NumberCharged-1)*(NumberParticles-NumberCharged)*(NumberParticles-NumberCharged-1))/static_cast<G4double>(denominator);  //re-recorrected JMQ 03/10/07

    return rj;
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
 

