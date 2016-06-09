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
// by V. Lara
//
//J. M. Quesada (July 08) 


#ifndef G4PreCompoundAlpha_h
#define G4PreCompoundAlpha_h 1

#include "G4PreCompoundIon.hh"
#include "G4ReactionProduct.hh"
#include "G4Alpha.hh"
#include "G4AlphaCoulombBarrier.hh"
#include "G4PreCompoundParameters.hh"

class G4PreCompoundAlpha : public G4PreCompoundIon
{
public:
  // default constructor
  G4PreCompoundAlpha():G4PreCompoundIon(4,2,&theAlphaCoulombBarrier,"Alpha") {}

  // copy constructor
  G4PreCompoundAlpha(const G4PreCompoundAlpha &right): G4PreCompoundIon(right){}

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


  G4ReactionProduct * GetReactionProduct() const;

private:

  virtual G4double GetRj(const G4int NumberParticles, const G4int NumberCharged);

  virtual G4double CrossSection(const  G4double K) ; 

  virtual G4double FactorialFactor(const G4double N, const G4double P);

  virtual G4double CoalescenceFactor(const G4double A);

  G4double GetOpt0(const G4double K);
  G4double GetOpt12(const G4double K);
  G4double GetOpt34(const G4double K);

  G4double GetAlpha();
  
  G4double GetBeta();

//data members

      G4AlphaCoulombBarrier theAlphaCoulombBarrier;
      G4double ResidualA;
      G4double ResidualZ; 
      G4double theA;
      G4double theZ;
      G4double ResidualAthrd;
      G4double FragmentA;
      G4double FragmentAthrd;

};
#endif





 



