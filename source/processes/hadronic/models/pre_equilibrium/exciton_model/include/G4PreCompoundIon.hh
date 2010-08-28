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
// $Id: G4PreCompoundIon.hh,v 1.8 2010-08-28 15:16:55 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// J. M. Quesada (August 2008).  
// Based  on previous work by V. Lara
//
// Modified:
// 20.08.2010 V.Ivanchenko added int Z and A and cleanup; added 
//                        G4ParticleDefinition to constructor,
//                        moved constructor and destructor to source,
//                        added inline methods

#ifndef G4PreCompoundIon_h
#define G4PreCompoundIon_h 1

#include "G4PreCompoundFragment.hh"

class G4PreCompoundIon : public G4PreCompoundFragment
{
public:

  G4PreCompoundIon(const G4ParticleDefinition*,
		   G4VCoulombBarrier * aCoulombBarrier);
  
  virtual ~G4PreCompoundIon();

protected:
          
  virtual G4double 
  ProbabilityDistributionFunction(G4double eKin, 
				  const G4Fragment& aFragment);

  virtual G4double CrossSection(G4double ekin) = 0; 

  virtual G4double 
  GetRj(G4int NumberParticles, G4int NumberCharged) = 0; 

  virtual G4double FactorialFactor(G4int N, G4int P) = 0;

  virtual G4double CoalescenceFactor(G4int A) = 0; 

  virtual G4double GetAlpha() = 0;

  inline G4double GetBeta();

  inline G4double GetOpt0(G4double ekin);

private:

  // default constructor
  G4PreCompoundIon();
  // operators
  G4PreCompoundIon(const G4PreCompoundIon &right);
  const G4PreCompoundIon& 
  operator= (const G4PreCompoundIon &right);
  G4int operator==(const G4PreCompoundIon &right) const;
  G4int operator!=(const G4PreCompoundIon &right) const;    

  G4double fact;
};

inline G4double G4PreCompoundIon::GetBeta()
{
  return -GetCoulombBarrier();
}

// *********************** OPT=0 : Dostrovski's cross section  ***************
inline G4double G4PreCompoundIon::GetOpt0(G4double K)
{
  G4double r0 = theParameters->Getr0()*ResidualA13();
  // cross section is now given in mb (r0 is in mm) for the sake of consistency
  //with the rest of the options
  return 1.e+25*CLHEP::pi*r0*r0*ResidualA13()*GetAlpha()*(1.+GetBeta()/K);
}

#endif
