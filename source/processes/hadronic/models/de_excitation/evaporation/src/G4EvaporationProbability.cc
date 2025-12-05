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
#include "G4InterfaceToXS.hh"
#include "G4IsotopeList.hh"
#include "G4DeexPrecoUtility.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

namespace
{
  const G4double explim = 160.;        // limit of G4Exp argument, OPTxs = 0
  const G4double lim = 2*CLHEP::MeV;   // limit on x-section interface, OPTxs = 1
  const G4double kmin = 20*CLHEP::keV; // low-energy limit on primary kinetic energy
}

G4EvaporationProbability::G4EvaporationProbability(G4int anA, G4int aZ, 
						   G4double aGamma) 
  : G4VEmissionProbability(aZ, anA), fGamma(aGamma)
{
  resA13 = lastA = muu = freeU = a0 = a1 = delta0 = delta1 = 0.0;
  pcoeff = fGamma*pEvapMass*CLHEP::millibarn
    /((CLHEP::pi*CLHEP::hbarc)*(CLHEP::pi*CLHEP::hbarc)); 

  if (1 == theZ && 1 == theA) { index = 1; }
  else if (1 == theZ && 2 == theA) { index = 2; }
  else if (1 == theZ && 3 == theA) { index = 3; }
  else if (2 == theZ && 3 == theA) { index = 4; }
  else if (2 == theZ && 4 == theA) { index = 5; }

  if (OPTxs == 1) {
    const G4ParticleDefinition* part = nullptr;
    if (index == 1) { part = G4Proton::Proton(); }
    else if (index == 2) { part = G4Deuteron::Deuteron(); }
    else if (index == 3) { part = G4Triton::Triton(); }
    else if (index == 4) { part = G4He3::He3(); }
    else if (index == 5) { part = G4Alpha::Alpha(); }
    else { part = G4Neutron::Neutron(); }
    fXSection = new G4InterfaceToXS(part, index);
  }
  
  if (0 == aZ) {
    ResetIntegrator(0.15*CLHEP::MeV, 0.01);
  } else {
    ResetIntegrator(0.20*CLHEP::MeV, 0.01);
  }
}

G4EvaporationProbability::~G4EvaporationProbability()
{
  delete fXSection;
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
  freeU = exEnergy;
  a0 = pNuclearLevelData->GetLevelDensity(fragZ, fragA, freeU);
  delta0 = pNuclearLevelData->GetPairingCorrection(fragZ, fragA);
  delta1 = pNuclearLevelData->GetPairingCorrection(resZ, resA);
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
    a1 = pNuclearLevelData->GetLevelDensity(resZ,resA,freeU);
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

G4double G4EvaporationProbability::ComputeProbability(G4double kinE, G4double CB)
{
  G4double K = std::max(kinE, kmin);
  // abnormal case - should never happens
  if(pMass < pEvapMass + pResMass + K) { return 0.0; }
    
  G4double K1 = pMass - pEvapMass - K;
  G4double mres = std::sqrt(K1*K1 - K*(2*pEvapMass + K));

  G4double excRes = mres - pResMass;
  if (excRes < 0.0) { return 0.0; }
  G4double K2 = 0.5*(pMass + pEvapMass + mres)*(pMass - pEvapMass - mres)/mres;
  G4double xs = CrossSection(K2, CB);
  if (xs <= 0.0) { return 0.0; }

  a1 = pNuclearLevelData->GetLevelDensity(resZ, resA, excRes);
  G4double E0 = std::max(freeU - delta0, 0.0);
  G4double E1 = std::max(excRes - delta1, 0.0);
  G4double prob = pcoeff*G4Exp(2.0*(std::sqrt(a1*E1) - std::sqrt(a0*E0)))*K*xs;
  return prob;
}

G4double 
G4EvaporationProbability::CrossSection(G4double kine, G4double CB)
{
  G4double K = std::max(kine, kmin);
  // compute power once
  if (OPTxs > 1 && 0 < index && resA != lastA) {
    lastA = resA;
    muu = G4KalbachCrossSection::ComputePowerParameter(resA, index);
  }
  // In the case of OPTxs = 0 this method is not called
  if (OPTxs == 1) {
    G4int Z = std::min(resZ, ZMAXNUCLEARDATA);
    if (0 == index) {
      G4double e1 = lowEnergyLimitMeV[Z];
      if (e1 == 0.0) { e1 = lim; }
      K = std::max(K, e1);
    } else {
      if (K < 0.5*CB) {
	recentXS = 0.0;
	return recentXS;
      }
      K = std::max(K, 2*CB);
    }
    G4double corr = G4DeexPrecoUtility::CorrectionFactor(index, theZ, resA13, CB, kine);
    recentXS = corr*fXSection->GetElementCrossSection(K, Z)/CLHEP::millibarn;

  } else if (OPTxs == 2) { 
    recentXS = G4ChatterjeeCrossSection::ComputeCrossSection(K, CB, resA13, muu,
                                                             index, theZ, resA); 
  } else if (OPTxs == 3) {
    recentXS = G4KalbachCrossSection::ComputeCrossSection(K, CB, resA13, muu,
                                                          index, theZ, theA, resA);
  }
  return recentXS;
}
