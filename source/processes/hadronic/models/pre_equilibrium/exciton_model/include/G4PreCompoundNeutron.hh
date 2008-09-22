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
//J. M. Quesada (July  08) 


#ifndef G4PreCompoundNeutron_h
#define G4PreCompoundNeutron_h 1

#include "G4PreCompoundNucleon.hh"
#include "G4ReactionProduct.hh"
#include "G4Neutron.hh"
#include "G4PreCompoundParameters.hh"
#include "Randomize.hh"
#include "G4NeutronCoulombBarrier.hh"


class G4PreCompoundNeutron : public G4PreCompoundNucleon
{
public:
  // default constructor
  G4PreCompoundNeutron() : G4PreCompoundNucleon(1,0,&theNeutronCoulombBarrier,"Neutron") {}

  // copy constructor
  G4PreCompoundNeutron(const G4PreCompoundNeutron &right): G4PreCompoundNucleon(right) {}

  // destructor
  ~G4PreCompoundNeutron() {}

  // operators  
  const G4PreCompoundNeutron & operator=(const G4PreCompoundNeutron &right) {
    if (&right != this) this->G4PreCompoundNucleon::operator=(right);
    return *this;
  }

  G4bool operator==(const G4PreCompoundNeutron &right) const
  { return G4PreCompoundNucleon::operator==(right);}
  
  G4bool operator!=(const G4PreCompoundNeutron &right) const
  { return G4PreCompoundNucleon::operator!=(right);}


  G4ReactionProduct * GetReactionProduct() const;
  
private:

  virtual G4double GetRj(const G4int NumberParticles, const G4int NumberCharged); 

  virtual G4double CrossSection(const  G4double K); 


  G4double GetOpt0(const G4double K);
  G4double GetOpt12(const G4double K);
  G4double GetOpt34(const G4double K);

  
  G4double GetAlpha();
  
  G4double GetBeta();
  
//data members

      G4NeutronCoulombBarrier theNeutronCoulombBarrier;
      G4double ResidualA;
      G4double ResidualZ; 
      G4double theA;
      G4double theZ;
      G4double ResidualAthrd;
      G4double FragmentA;
      G4double FragmentAthrd;

 
};

 
#endif

