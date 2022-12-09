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
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICB is set true, by default is false) 
//
// JMQ (14 february 2009) bug fixed in emission width: hbarc instead of 
//                        hbar_Planck in the denominator
//
// V.Ivanchenko general clean-up since 2010
//
#include "G4EvaporationProbability.hh"
#include "G4NuclearLevelData.hh"
#include "G4VCoulombBarrier.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PairingCorrection.hh"
#include "G4NucleiProperties.hh"
#include "G4KalbachCrossSection.hh"
#include "G4ChatterjeeCrossSection.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

static const G4double explim = 160.;

G4EvaporationProbability::G4EvaporationProbability(G4int anA, G4int aZ, 
						   G4double aGamma) 
  : G4VEmissionProbability(aZ, anA), fGamma(aGamma)
{
  resA13 = lastA = muu = freeU = a0 = delta1 = 0.0;
  pcoeff = fGamma*pEvapMass*CLHEP::millibarn
    /((CLHEP::pi*CLHEP::hbarc)*(CLHEP::pi*CLHEP::hbarc)); 

  if(0 == theZ)      { index = 0; }
  else if(1 == theZ) { index = theA; }
  else               { index = theA + 1; }
  if(0 == aZ) {
    ResetIntegrator(30, 0.25*CLHEP::MeV, 0.02);
  } else {
    ResetIntegrator(30, 0.5*CLHEP::MeV, 0.03);
  }
}

G4double G4EvaporationProbability::CalcAlphaParam(const G4Fragment&)
{
  return 1.0;
}
 
G4double G4EvaporationProbability::CalcBetaParam(const G4Fragment&)
{
  return 1.0;
}

G4double G4EvaporationProbability::TotalProbability(
  const G4Fragment& fragment, G4double minEnergy, G4double maxEnergy, 
  G4double CB, G4double exEnergy)
{
  G4int fragA = fragment.GetA_asInt();
  G4int fragZ = fragment.GetZ_asInt();
  G4double U = fragment.GetExcitationEnergy(); 
  a0 = pNuclearLevelData->GetLevelDensity(fragZ,fragA,U);
  freeU = exEnergy;
  delta1 = pNuclearLevelData->GetPairingCorrection(resZ,resA);
  resA13 = pG4pow->Z13(resA);
  /*      
  G4cout << "G4EvaporationProbability: Z= " << theZ << " A= " << theA 
	 << " resZ= " << resZ << " resA= " << resA 
	 << " fragZ= " << fragZ << " fragA= " << fragA 
	 << "\n   freeU= " << freeU  
	 << " a0= " << a0 << " OPT= " << OPTxs << " emin= " 
         << minEnergy << " emax= " << maxEnergy 
	 << " CB= " << CB << G4endl;
  */
  if (OPTxs==0) {

    G4double SystemEntropy = 2.0*std::sqrt(a0*freeU);
    const G4double RN2 = 2.25*CLHEP::fermi*CLHEP::fermi
      /(CLHEP::twopi*CLHEP::hbar_Planck*hbar_Planck);

    G4double Alpha = CalcAlphaParam(fragment);
    G4double Beta = CalcBetaParam(fragment);

    // to be checked where to use a0, where - a1	
    G4double a1 = pNuclearLevelData->GetLevelDensity(resZ,resA,freeU);
    G4double GlobalFactor = fGamma*Alpha*pEvapMass*RN2*resA13*resA13/(a1*a1);
    
    G4double maxea = maxEnergy*a1;
    G4double Term1 = Beta*a1 - 1.5 + maxea;
    G4double Term2 = (2.0*Beta*a1-3.0)*std::sqrt(maxea) + 2*maxea;
	
    G4double ExpTerm1 = (SystemEntropy <= explim) ? G4Exp(-SystemEntropy) : 0.0;
	
    G4double ExpTerm2 = 2.*std::sqrt(maxea) - SystemEntropy;
    ExpTerm2 = std::min(ExpTerm2, explim);
    ExpTerm2 = G4Exp(ExpTerm2);
	
    pProbability = GlobalFactor*(Term1*ExpTerm1 + Term2*ExpTerm2);
             
  } else {
    // if Coulomb barrier cutoff is superimposed for all cross sections 
    // then the limit is the Coulomb Barrier
    pProbability = IntegrateProbability(minEnergy, maxEnergy, CB);
  }
  /*
  G4cout << "TotalProbability:  Emin=" << minEnergy << " Emax= " << maxEnergy 
	 << " CB= " << CB << " prob=" << pProbability << G4endl;
  */
  return pProbability;
}

G4double G4EvaporationProbability::ComputeProbability(G4double K, G4double CB)
{ 
  G4double E0 = freeU;
  // abnormal case - should never happens
  if(pMass < pEvapMass + pResMass) { return 0.0; }
    
  G4double m02 = pMass*pMass;
  G4double m12 = pEvapMass*pEvapMass;
  G4double mres = std::sqrt(m02 + m12 - 2.*pMass*(pEvapMass + K));

  G4double excRes = mres - pResMass;
  G4double E1 = excRes - delta1;
  if(E1 <= 0.0) { return 0.0; }
  G4double a1 = pNuclearLevelData->GetLevelDensity(resZ,resA,excRes);
  G4double xs = CrossSection(K, CB); 
  G4double prob = pcoeff*G4Exp(2.0*(std::sqrt(a1*E1) - std::sqrt(a0*E0)))*K*xs;
  return prob;
}

G4double 
G4EvaporationProbability::CrossSection(G4double K, G4double CB)
{
  // compute power once
  if(resA != lastA) {
    lastA = resA;
    if(0 < index) 
      muu = G4KalbachCrossSection::ComputePowerParameter(resA, index);
  }
  G4double res = 0.0;
  if(OPTxs <= 2) { 
    res = G4ChatterjeeCrossSection::ComputeCrossSection(K, CB, resA13, muu,
							index, theZ, resA); 
  } else {
    // added barrier penetration factor
    G4double elim = 0.6*CB;
    if(K > elim) {
      res = G4KalbachCrossSection::ComputeCrossSection(K, elim, resA13, muu,
						       index, theZ, theA, resA);
      res *= (1.0 - elim/K);
    }
  }
  return res;
}
