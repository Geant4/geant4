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
// $Id: G4PreCompoundDeuteron.hh,v 1.11 2008-07-24 13:53:32 quesada Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara
//
//
//J. M. Quesada (Dic. 07) Added combinatorial factor Rj. Removed  GetAlpha and GetBeta methods

#ifndef G4PreCompoundDeuteron_h
#define G4PreCompoundDeuteron_h 1

#include "G4PreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Deuteron.hh"
#include "G4DeuteronCoulombBarrier.hh"



class G4PreCompoundDeuteron : public G4PreCompoundIon
{
public:
  // default constructor
  G4PreCompoundDeuteron():G4PreCompoundIon(2,1,&theDeuteronCoulombBarrier,"Deuteron") {}

  // copy constructor
  G4PreCompoundDeuteron(const G4PreCompoundDeuteron &right): G4PreCompoundIon(right) {}

  // destructor
  ~G4PreCompoundDeuteron() {}

  // operators  
  const G4PreCompoundDeuteron & operator=(const G4PreCompoundDeuteron &right) {
    if (&right != this) this->G4PreCompoundIon::operator=(right);
    return *this;
  }

  G4bool operator==(const G4PreCompoundDeuteron &right) const
  { return G4PreCompoundIon::operator==(right);}

  
  G4bool operator!=(const G4PreCompoundDeuteron &right) const
  { return G4PreCompoundIon::operator!=(right);}


  G4ReactionProduct * GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct =
      new G4ReactionProduct(G4Deuteron::DeuteronDefinition());
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
    G4double denominator = NumberParticles*(NumberParticles-1);
    if(NumberCharged >= 1 && (NumberParticles-NumberCharged) >= 1) rj = 2.0*static_cast<G4double>(NumberCharged*(NumberParticles-NumberCharged))/static_cast<G4double>(denominator); //JMQ re-recorrected 03/10/07
 

    return rj;
  }




  virtual G4double FactorialFactor(const G4double N, const G4double P)
  {
    return (N-1.0)*(N-2.0)*(P-1.0)*P/2.0;
  }
  
  virtual G4double CoalescenceFactor(const G4double A)
  {
    return 16.0/A;
  }    
private:

  G4DeuteronCoulombBarrier theDeuteronCoulombBarrier;

};

#endif
 

