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

// $Id: G4PreCompoundProton.hh,v 1.11 2008-07-24 13:53:32 quesada Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara
//
//J. M. Quesada (Dic. 07) Added combinatorial factor Rj. Removed  GetAlpha and GetBeta methods

#ifndef G4PreCompoundProton_h
#define G4PreCompoundProton_h 1

#include "G4PreCompoundNucleon.hh"
#include "G4ReactionProduct.hh"
#include "G4Proton.hh"
#include "G4PreCompoundParameters.hh"
#include "Randomize.hh"
#include "G4ProtonCoulombBarrier.hh"

class G4PreCompoundProton : public G4PreCompoundNucleon
{
public:
  // default constructor
  G4PreCompoundProton():G4PreCompoundNucleon(1,1,&theProtonCoulombBarrier,"Proton") {}

  // copy constructor
  G4PreCompoundProton(const G4PreCompoundProton &right): G4PreCompoundNucleon(right) {}

  // destructor
  ~G4PreCompoundProton() {};

  // operators  
  const G4PreCompoundProton & operator=(const G4PreCompoundProton &right) {
    if (&right != this) this->G4PreCompoundNucleon::operator=(right);
    return *this;
  };

  G4bool operator==(const G4PreCompoundProton &right) const
  { return G4PreCompoundNucleon::operator==(right);}

  
  G4bool operator!=(const G4PreCompoundProton &right) const
  { return G4PreCompoundNucleon::operator!=(right);}


  G4ReactionProduct * GetReactionProduct() const
  {
    G4ReactionProduct * theReactionProduct = 
      new G4ReactionProduct(G4Proton::ProtonDefinition());
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
    if(NumberParticles > 0) rj = static_cast<G4double>(NumberCharged)/static_cast<G4double>(NumberParticles);
    return rj;
  }

  
  virtual G4bool IsItPossible(const G4Fragment& aFragment)
  {
    return (aFragment.GetNumberOfCharged() >= 1);
  }
private:

  G4ProtonCoulombBarrier theProtonCoulombBarrier;

  
};

#endif
 
