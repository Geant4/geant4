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
// $Id: G4GEMProbability.hh 103162 2017-03-20 09:40:58Z gcosmo $
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

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4VCoulombBarrier.hh"
#include "G4PairingCorrection.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"

class G4GEMProbability : public G4VEmissionProbability
{
public:

  G4GEMProbability(G4int anA, G4int aZ, G4double aSpin);
    
  virtual ~G4GEMProbability();

  // not used for evaporation
  virtual G4double EmissionProbability(const G4Fragment& fragment,
				       G4double maxKineticEnergy);

  void Dump() const;

  inline G4double GetSpin(void) const;

  inline void SetCoulomBarrier(const G4VCoulombBarrier * aCoulombBarrierStrategy);

  inline G4double GetCoulombBarrier(const G4Fragment& fragment) const; 

  inline G4double CalcAlphaParam(const G4Fragment & ) const;

  inline G4double CalcBetaParam(const G4Fragment & ) const;
        
private:
    
  G4double CalcProbability(const G4Fragment & fragment, 
			   G4double MaximalKineticEnergy,
			   G4double V);

  inline G4double CCoeficient(G4int) const;

  inline G4double I0(G4double t);
  inline G4double I1(G4double t, G4double tx);
  inline G4double I2(G4double s0, G4double sx);
  G4double I3(G4double s0, G4double sx);

  // Copy constructor
  G4GEMProbability();
  G4GEMProbability(const G4GEMProbability &right);    
  const G4GEMProbability & operator=(const G4GEMProbability &right);
  G4bool operator==(const G4GEMProbability &right) const;
  G4bool operator!=(const G4GEMProbability &right) const;
    
  // Data Members
  G4Pow*   fG4pow;
  G4PairingCorrection* fPairCorr;
    
  G4VLevelDensityParameter * theEvapLDPptr;
    
  // Spin is fragment spin
  G4double Spin;

  // Coulomb Barrier
  const G4VCoulombBarrier * theCoulombBarrierPtr;
  
protected:

  G4double fPlanck;

  // Resonances Energy
  std::vector<G4double> ExcitEnergies;
    
  // Resonances Spin 
  std::vector<G4double> ExcitSpins;

  // Resonances half lifetime
  std::vector<G4double> ExcitLifetimes;

};

inline G4double G4GEMProbability::GetSpin(void) const 
{ 
  return Spin; 
}

inline void 
G4GEMProbability::SetCoulomBarrier(const G4VCoulombBarrier * aCoulombBarrierStrategy)
{
  theCoulombBarrierPtr = aCoulombBarrierStrategy;
}

inline G4double 
G4GEMProbability::GetCoulombBarrier(const G4Fragment& fragment) const 
{
  G4double res = 0.0;
  if (theCoulombBarrierPtr) {
    G4int Acomp = fragment.GetA_asInt();
    G4int Zcomp = fragment.GetZ_asInt();
    res = theCoulombBarrierPtr->GetCoulombBarrier(Acomp-theA, Zcomp-theZ,
				fragment.GetExcitationEnergy() -
				fPairCorr->GetPairingCorrection(Acomp,Zcomp));
  }
  return res;
}

inline G4double G4GEMProbability::CCoeficient(G4int aZ) const
{
  //JMQ 190709 C's values from Furihata's paper 
  //(notes added on proof in Dostrovskii's paper) 
  //data = {{20, 0.}, {30, -0.06}, {40, -0.10}, {50, -0.10}};
  G4double C = 0.0;
  if (aZ >= 50){
    C=-0.10/G4double(theA);
  } else if (aZ > 20) {
    C=(0.123482-0.00534691*aZ-0.0000610624*aZ*aZ+5.93719*1e-7*aZ*aZ*aZ+
       1.95687*1e-8*aZ*aZ*aZ*aZ)/G4double(theA);
  }
  return C;
}


inline G4double G4GEMProbability::CalcAlphaParam(const G4Fragment & fragment) const
{
  //JMQ 190709 values according to Furihata's paper (based on notes added 
  //on proof in Dostrovskii's paper)
  G4double res;
  if(theZ == 0) {
    res = 0.76+1.93/fG4pow->Z13(fragment.GetA_asInt()-theA);
  } else {
    res = 1.0 + CCoeficient(fragment.GetZ_asInt()-theZ);
  }
  return res;
}

inline G4double 
G4GEMProbability::CalcBetaParam(const G4Fragment & fragment) const
{
  //JMQ 190709 values according to Furihata's paper (based on notes added 
  //on proof in Dostrovskii's paper)
  G4double res;
  if(theZ == 0) {
    res = (1.66/fG4pow->Z23(fragment.GetA_asInt()-theA)-0.05)*CLHEP::MeV/
	   CalcAlphaParam(fragment);
  } else {
    res = -GetCoulombBarrier(fragment);
  }
  return res;
}

inline G4double G4GEMProbability::I0(G4double t)
{
  return G4Exp(t) - 1.0;
}

inline G4double G4GEMProbability::I1(G4double t, G4double tx)
{
  return (t - tx + 1.0)*G4Exp(tx) - t - 1.0;
}


inline G4double G4GEMProbability::I2(G4double s0, G4double sx)
{
  G4double S = 1.0/std::sqrt(s0);
  G4double Sx = 1.0/std::sqrt(sx);
  
  G4double p1 = S*S*S*( 1.0 + S*S*( 1.5 + 3.75*S*S) );
  G4double p2 = Sx*Sx*Sx*( 1.0 + Sx*Sx*( 1.5 + 3.75*Sx*Sx) )*G4Exp(sx-s0);
  
  return p1-p2;
}


#endif
