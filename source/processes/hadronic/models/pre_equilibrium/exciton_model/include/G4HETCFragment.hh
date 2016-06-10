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
// $Id: G4HETCFragment.hh 90337 2015-05-26 08:34:27Z gcosmo $
//
// by V. Lara
//
// Modified:  
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup

#ifndef G4HETCFragment_h
#define G4HETCFragment_h 1

#include "G4VPreCompoundFragment.hh"
#include "Randomize.hh"

class G4HETCFragment : public G4VPreCompoundFragment
{
public:  

  G4HETCFragment(const G4ParticleDefinition*,
		 G4VCoulombBarrier * aCoulombBarrier);
  
  virtual ~G4HETCFragment();

  G4double CalcEmissionProbability(const G4Fragment & aFragment);

protected:

  virtual G4double K(const G4Fragment & aFragment) = 0;
    
  virtual G4double GetSpinFactor() const = 0;
  virtual G4double GetAlpha() const = 0;
  virtual G4double GetBeta() const = 0;

  inline G4double BetaRand(G4int N, G4int L) const;
  
private:

  // This method performs integration for probability function over 
  // fragment kinetic energy
  G4double IntegrateEmissionProbability(G4double & Low, 
					G4double & Up, 
					const G4Fragment & aFragment);	

  G4HETCFragment();
  G4HETCFragment(const G4HETCFragment &right);
  const G4HETCFragment& 
  operator= (const G4HETCFragment &right);  
  G4int operator==(const G4HETCFragment &right) const;
  G4int operator!=(const G4HETCFragment &right) const;

  G4double r2norm;
};

inline G4double G4HETCFragment::BetaRand(G4int N, G4int L) const
{
  G4double Y1 = G4RandGamma::shoot(N,1);
  G4double Y2 = G4RandGamma::shoot(L,1);
  
  return Y1/(Y1+Y2);
}

#endif
