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
// $Id: G4PreCompoundTriton.hh,v 1.10 2008-07-24 13:53:32 quesada Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara
//
//J. M. Quesada (Dic. 07) Added combinatorial factor Rj. Removed  GetAlpha and GetBeta methods

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

//JMQ (Sep. 07) combinatorial factor Rj
  virtual G4double GetRj(const G4int NumberParticles, const G4int NumberCharged)
  {
    G4double rj = 0.0;
    G4double denominator = NumberParticles*(NumberParticles-1)*(NumberParticles-2);
    if(NumberCharged >= 1 && (NumberParticles-NumberCharged) >= 2) rj = 3.0*static_cast<G4double>(NumberCharged*(NumberParticles-NumberCharged)*(NumberParticles-NumberCharged-1))/static_cast<G4double>(denominator); //re-recorrected JMQ 03/10/07

    return rj;
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
