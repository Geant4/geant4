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
//J.M. Quesada (August2008). Based on:
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

using namespace std;

G4EvaporationProbability::G4EvaporationProbability(G4int anA, G4int aZ, 
						   G4double aGamma, 
						   G4VCoulombBarrier*) 
  : G4VEmissionProbability(aZ, anA), Gamma(aGamma)
{
  resZ = resA = fragA = fragZ = 0;
  resA13 = muu = resMass = Mass = U = delta0 = delta1 = a0 = probmax = 0.0;
  partMass = G4NucleiProperties::GetNuclearMass(theA, theZ);
  pcoeff = Gamma*partMass*CLHEP::millibarn
    /((CLHEP::pi*CLHEP::hbarc)*(CLHEP::pi*CLHEP::hbarc)); 

  if(0 == theZ)      { index = 0; }
  else if(1 == theZ) { index = theA; }
  else               { index = theA + 1; }
  fLevelData = G4NuclearLevelData::GetInstance(); 
  if(0 == aZ) {
    ResetIntegrator(30, 0.25*CLHEP::MeV, 0.02);
  } else if(1 == aZ && 1 == anA) {
    ResetIntegrator(20, 0.5*CLHEP::MeV, 0.03);
  } else {
    ResetIntegrator(20, CLHEP::MeV, 0.04);
  }
}

G4EvaporationProbability::~G4EvaporationProbability() 
{}

G4double G4EvaporationProbability::EmissionProbability(
         const G4Fragment&, G4double)
{
  return 0.0;
}

G4double G4EvaporationProbability::TotalProbability(
  const G4Fragment & fragment, G4double minEnergy, G4double maxEnergy, 
  G4double CoulombBarrier)
{
  fragA = fragment.GetA_asInt();
  fragZ = fragment.GetZ_asInt();
  resA  = fragA - theA;
  resZ  = fragZ - theZ;

  G4double fragMass = fragment.GetGroundStateMass();
  U = fragment.GetExcitationEnergy();
  Mass = fragMass + U;
  delta0 = std::max(0.0, fPairCorr->GetPairingCorrection(fragA,fragZ));
  delta1 = std::max(0.0, fPairCorr->GetPairingCorrection(resA,resZ));
  resMass = G4NucleiProperties::GetNuclearMass(resA, resZ);
  resA13 = fG4pow->Z13(resA);
  a0 = LevelDensity*fragA;
  /*      
  G4cout << "G4EvaporationProbability: Z= " << theZ << " A= " << theA 
	 << " resZ= " << resZ << " resA= " << resA 
	 << " fragZ= " << fragZ << " fragA= " << fragA 
	 << "\n      U= " << U << " d0= " << delta0 << " a0= " << a0 
	 << " OPT= " << OPTxs << G4endl;
  */
  if(U < delta0 || maxEnergy <= minEnergy) { return 0.0; }
   
  G4double Width = 0.0;
  if (OPTxs==0) {

    G4double SystemEntropy = 2.0*std::sqrt(a0*(U-delta0));
								  
    static const G4double RN2 = 
      2.25*fermi*fermi/(twopi* hbar_Planck*hbar_Planck);

    G4double Alpha = CalcAlphaParam(fragment);
    G4double Beta = CalcBetaParam(fragment);

    // to be checked where to use a0, where - a1	
    G4double a1 = LevelDensity*resA;
    G4double GlobalFactor = Gamma*Alpha*partMass*RN2*resA13*resA13/(a1*a1);
    
    G4double maxea = maxEnergy*a1;
    G4double Term1 = Beta*a1 - 1.5 + maxea;
    G4double Term2 = (2.0*Beta*a1-3.0)*std::sqrt(maxea) + 2*maxea;
	
    static const G4double explim = 350.;
    G4double ExpTerm1 = (SystemEntropy <= explim) ? G4Exp(-SystemEntropy) : 0.0;
	
    G4double ExpTerm2 = 2.*std::sqrt(maxea) - SystemEntropy;
    ExpTerm2 = std::min(ExpTerm2, explim);
    ExpTerm2 = G4Exp(ExpTerm2);
	
    Width = GlobalFactor*(Term1*ExpTerm1 + Term2*ExpTerm2);
             
  } else {

    // compute power once
    if(OPTxs <= 2) { 
      muu =  G4ChatterjeeCrossSection::ComputePowerParameter(resA, index);
    } else {
      muu = G4KalbachCrossSection::ComputePowerParameter(resA, index);
    }
    // if Coulomb barrier cutoff is superimposed for all cross sections 
    // then the limit is the Coulomb Barrier
    //Width = IntegrateEmissionProbability(minEnergy, maxEnergy, 
    //					 CoulombBarrier);
    Width = IntegrateProbability(minEnergy, maxEnergy, CoulombBarrier);
  }
  return Width;
}

G4double G4EvaporationProbability::ComputeProbability(G4double K, G4double cb)
{ 
  //G4cout << "### G4EvaporationProbability::ProbabilityDistributionFunction" 
  //  << G4endl;

  G4double E0 = U - delta0;
  //G4double E1 = Mass - partMass - resMass - delta1 - K;
  G4double E1 = sqrt((Mass - partMass)*(Mass - partMass) - 2*Mass*K) 
    - resMass - delta1;
  /*  
  G4cout << "PDF: FragZ= " << fragZ << " FragA= " << fragA
  	 << " Z= " << theZ << "  A= " << theA 
	 << " K= " << K << " E0= " << E0 << " E1= " << E1 << G4endl;
  */
  if(E1 < 0.0) { return 0.0; }

  G4double a1 = LevelDensity*resA;
  G4double Prob = pcoeff*G4Exp(2.0*(std::sqrt(a1*E1) - std::sqrt(a0*E0)))
    *K*CrossSection(K, cb);

  //G4cout << "Evap prob: " << Prob << "  fVerbose= " << fVerbose << G4endl;
  return Prob;
}

G4double 
G4EvaporationProbability::CrossSection(G4double K, G4double cb)
{
  G4double res;
  if(OPTxs <= 2) { 
    res = G4ChatterjeeCrossSection::ComputeCrossSection(K, cb, resA13, muu,
							index, theZ, resA); 
  } else { 
    res = G4KalbachCrossSection::ComputeCrossSection(K, cb, resA13, muu,
						     index, theZ, theA, resA);
  }
  //G4cout << "     K= " << K << "  res= " << res << " muu= " << muu << G4endl;
  return res;
}  

G4double 
G4EvaporationProbability::SampleKineticEnergy(G4double minKinEnergy, 
					      G4double maxKinEnergy,
					      G4double)
{
  /*
    G4cout << "### Sample probability Emin= " << minKinEnergy 
	   << " Emax= " << maxKinEnergy 
	   << "  Z= " << theZ << " A= " << theA << G4endl;
  */
  G4double T = 0.0;
  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  static const G4int nmax = 100;
  if (OPTxs==0) {
    // JMQ:
    // It uses Dostrovsky's approximation for the inverse reaction cross
    // in the probability for fragment emission
    // MaximalKineticEnergy energy in the original version (V.Lara) was 
    // calculated at the Coulomb barrier.
    
    G4double Rb = 4.0*a0*maxKinEnergy;
    G4double RbSqrt = std::sqrt(Rb);
    G4double PEX1 = 0.0;
    if (RbSqrt < 160.0) { PEX1 = G4Exp(-RbSqrt); }
    G4double Rk = 0.0;
    G4double FRk = 0.0;
    G4int nn = 0;
    do {
      G4double RandNumber = rndm->flat();
      Rk = 1.0 + (1./RbSqrt)*G4Log(RandNumber + (1.0-RandNumber)*PEX1);
      G4double Q1 = 1.0;
      G4double Q2 = 1.0;
      if (theZ == 0) { // for emitted neutron
        G4double Beta = (2.12/(resA13*resA13) - 0.05)*MeV/(0.76 + 2.2/resA13);
        Q1 = 1.0 + Beta/maxKinEnergy;
        Q2 = Q1*std::sqrt(Q1);
      } 
      
      static const G4double ssqr3  = 1.5*std::sqrt(3.0);
      FRk = ssqr3 * Rk * (Q1 - Rk*Rk)/Q2;
      if(nn > nmax) { break; }
      ++nn;
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (FRk < rndm->flat());
    
    T = maxKinEnergy * (1.0-Rk*Rk) + minKinEnergy;

  } else { 

    if(fVerbose > 1) {
      G4cout << "###=== SampleEnergy: " << " Z= " << theZ
	     << " A= " << theA << " FragZ= " << fragZ 
	     << " FragA= " << fragA << G4endl; 
    }
    T = SampleEnergy();

  }
  //G4cout << "-- new Z= " << theZ << " A= " << theA << " ekin= " << T << G4endl; 
  return T;
  //return fLevelData->FindLevel(resZ, resA, resMass, Mass, partMass, T);
}
