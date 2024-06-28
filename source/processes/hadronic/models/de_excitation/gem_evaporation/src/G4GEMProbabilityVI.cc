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
// GEM de-excitation model
// by V. Ivanchenko (July 2019)
//
#include "G4GEMProbabilityVI.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4PairingCorrection.hh"
#include "G4NucleiProperties.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"


G4GEMProbabilityVI::G4GEMProbabilityVI(G4int anA, G4int aZ, const G4LevelManager* p) 
  : G4VEmissionProbability(aZ, anA), lManager(p)
{
  fragA  = fragZ = 0;
  resA13 = U = delta0 = delta1 = a0 = a1 = probmax = alphaP = betaP = 0.0;
  Umax = bCoulomb = 0.0;
  Gamma = 1.0;
  pcoeff = Gamma*pEvapMass*CLHEP::millibarn
    /((CLHEP::pi*CLHEP::hbarc)*(CLHEP::pi*CLHEP::hbarc)); 
  coeff = CLHEP::fermi*CLHEP::fermi/(CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc);

  isExcited = (!lManager || 0.0 == lManager->MaxLevelEnergy()) ? false : true;
  A13 = pG4pow->Z13(theA);

  if(0 == aZ) {
    ResetIntegrator(30, 0.25*CLHEP::MeV, 0.02);
  } else {
    ResetIntegrator(30, 0.5*CLHEP::MeV, 0.03);
  }
}

G4double G4GEMProbabilityVI::TotalProbability(
                             const G4Fragment& fragment,
			     const G4double tmin, const G4double tmax, 
			     const G4double CB, const G4double exEnergy,
			     const G4double exEvap)
{
  fragA = fragment.GetA_asInt();
  fragZ = fragment.GetZ_asInt();

  bCoulomb = CB;
  U = fragment.GetExcitationEnergy();
  delta0 = pNuclearLevelData->GetPairingCorrection(fragZ,fragA);
  delta1 = pNuclearLevelData->GetPairingCorrection(resZ,resA);
  Umax = pMass - pEvapMass - pResMass - CB;
  if(0.0 >= Umax) { return 0.0; }

  resA13 = pG4pow->Z13(resA);
  a0 = pNuclearLevelData->GetLevelDensity(fragZ,fragA,U);
  const G4double twoMass = pMass + pMass;
  const G4double evapMass2 = pEvapMass*pEvapMass;
  G4double ekinmax = 
     ((pMass-pResMass)*(pMass+pResMass) + evapMass2)/twoMass - pEvapMass;
  G4double ekinmin = 
      std::max((CB*(twoMass - CB) + evapMass2)/twoMass - pEvapMass,0.0);
  if(ekinmax <= ekinmin) { return 0.0; }
  pProbability = IntegrateProbability(ekinmin, ekinmax, CB);
  pProbability += tmax - tmin + exEnergy -exEvap;
  /*  
  G4cout << "G4GEMProbabilityVI: Z= " << theZ << " A= " << theA 
	 << " resZ= " << resZ << " resA= " << resA 
	 << " fragZ= " << fragZ << " fragA= " << fragA 
         << " prob= " << pProbability 
	 << "\n   U= " << U << " Umax= " << Umax << " d0= " << delta0 
         << " a0= " << a0 << G4endl;
  */
  return pProbability;
}

G4double G4GEMProbabilityVI::ComputeProbability(G4double ekin, G4double)
{ 
  // abnormal case - should never happens
  if(pMass < pEvapMass + pResMass) { return 0.0; }
    
  const G4double m02   = pMass*pMass;
  const G4double m12   = pEvapMass*pEvapMass;
  const G4double mres  = std::sqrt(m02 + m12 - 2.*pMass*(pEvapMass + ekin));

  G4double excRes = std::max(mres - pResMass, 0.0);
  a1 = pNuclearLevelData->GetLevelDensity(resZ,resA,excRes);
  G4double prob = 0.5; //CrossSection(0.0, excRes);

  //G4cout<<"### G4GEMProbabilityVI::ComputeProbability: Ekin(MeV)= "<<ekin 
  //<< " excRes(MeV)= " << excRes << " prob= " << prob << << G4endl;
  return prob;
}

G4double G4GEMProbabilityVI::SampleEnergy(
                             const G4double tmin, const G4double tmax, 
			     const G4double CB, const G4double exEnergy,
			     const G4double exEvap)
{
  G4double ekin = tmax - tmin - CB -exEnergy + exEvap;
  return ekin;
}
