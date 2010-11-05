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
// $Id: G4GEMProbability.hh,v 1.6 2010-11-05 14:42:52 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
//
// Geant4 header G4GEMProbability
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept 2001) 
//
// 18.05.2010 V.Ivanchenko trying to speedup the most slow method
//            by usage of G4Pow, integer Z and A; moved constructor, 
//            destructor and virtual functions to source
//

#ifndef G4GEMProbability_h
#define G4GEMProbability_h 1

#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4VCoulombBarrier.hh"
#include "G4PairingCorrection.hh"

class G4Pow;

class G4GEMProbability : public G4VEmissionProbability
{
public:

  // Default constructor - should not be used
  G4GEMProbability();

  // Only available constructor
  G4GEMProbability(G4int anA, G4int aZ, G4double aSpin);
    
  virtual ~G4GEMProbability();

  inline G4int GetZ_asInt(void) const { return theZ; }
	
  inline G4int GetA_asInt(void) const { return theA;}
	
  inline G4double GetZ(void) const { return theZ; }
	
  inline G4double GetA(void) const { return theA;}

  inline G4double GetSpin(void) const { return Spin; }

  inline G4double GetNormalization(void) const { return Normalization; }
    
  inline void SetCoulomBarrier(const G4VCoulombBarrier * aCoulombBarrierStrategy)
  {
    theCoulombBarrierPtr = aCoulombBarrierStrategy;
  }

  inline G4double GetCoulombBarrier(const G4Fragment& fragment) const 
  {
    G4double res = 0.0;
    if (theCoulombBarrierPtr) 
      {
	G4int Acomp = fragment.GetA_asInt();
	G4int Zcomp = fragment.GetZ_asInt();
	res = theCoulombBarrierPtr->GetCoulombBarrier(Acomp-theA, Zcomp-theZ,
	  fragment.GetExcitationEnergy()-fPairCorr->GetPairingCorrection(Acomp,Zcomp));
      }
    return res;
  }
    
  virtual G4double CalcAlphaParam(const G4Fragment & ) const;
  virtual G4double CalcBetaParam(const G4Fragment & ) const;
    
protected:
  
  inline void SetExcitationEnergiesPtr(std::vector<G4double> * anExcitationEnergiesPtr) 
  {
    ExcitationEnergies = anExcitationEnergiesPtr;
  }
  
  inline void SetExcitationSpinsPtr(std::vector<G4double> * anExcitationSpinsPtr)
  {
    ExcitationSpins = anExcitationSpinsPtr;
  }

  inline void SetExcitationLifetimesPtr(std::vector<G4double> * anExcitationLifetimesPtr)
  {
    ExcitationLifetimes = anExcitationLifetimesPtr;
  }

private:

  // Copy constructor
  G4GEMProbability(const G4GEMProbability &right);
    
  const G4GEMProbability & operator=(const G4GEMProbability &right);
  G4bool operator==(const G4GEMProbability &right) const;
  G4bool operator!=(const G4GEMProbability &right) const;
    
public:

  G4double EmissionProbability(const G4Fragment & fragment, G4double anEnergy);
  
private:

  G4double CalcProbability(const G4Fragment & fragment, G4double MaximalKineticEnergy,
			   G4double V);

  virtual G4double CCoeficient(G4double ) const;

  inline G4double I0(G4double t);
  inline G4double I1(G4double t, G4double tx);
  inline G4double I2(G4double s, G4double sx);
  G4double I3(G4double s, G4double sx);
    
  // Data Members

  G4Pow*   fG4pow;
  G4PairingCorrection* fPairCorr;
    
  G4VLevelDensityParameter * theEvapLDPptr;
	
  G4int theA;
  G4int theZ;
    
  // Spin is fragment spin
  G4double Spin;

  // Coulomb Barrier
  const G4VCoulombBarrier * theCoulombBarrierPtr;
  
  // Resonances Energy
  std::vector<G4double> * ExcitationEnergies;
    
  // Resonances Spin 
  std::vector<G4double> * ExcitationSpins;

  // Resonances half lifetime
  std::vector<G4double> * ExcitationLifetimes;

  // Normalization
  G4double Normalization;
    
};

inline G4double G4GEMProbability::I0(G4double t)
{
  return std::exp(t) - 1.0;
}

inline G4double G4GEMProbability::I1(G4double t, G4double tx)
{
  return (t - tx + 1.0)*std::exp(tx) - t - 1.0;
}


inline G4double G4GEMProbability::I2(G4double s, G4double sx)
{
  G4double S = 1.0/std::sqrt(s);
  G4double Sx = 1.0/std::sqrt(sx);
  
  G4double p1 = S*S*S*( 1.0 + S*S*( 1.5 + 3.75*S*S) );
  G4double p2 = Sx*Sx*Sx*( 1.0 + Sx*Sx*( 1.5 + 3.75*Sx*Sx) )*std::exp(sx-s);
  
  return p1-p2;
}


#endif
