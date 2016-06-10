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
// $Id: G4PhotonEvaporation.cc 94676 2015-12-02 09:51:20Z gunter $
//
// -------------------------------------------------------------------
//
//      GEANT4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PhotonEvaporation
//
//      Author:        Vladimir Ivantchenko
//
//      Creation date: 20 December 2011
//
//Modifications:
//
// 
// -------------------------------------------------------------------
//

#include "G4PhotonEvaporation.hh"

#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"
#include "G4FragmentVector.hh"
#include "G4GammaTransition.hh"
#include "G4PolarizedGammaTransition.hh"
#include "G4Pow.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

G4float G4PhotonEvaporation::GREnergy[] = {0.0f};
G4float G4PhotonEvaporation::GRWidth[] = {0.0f};
const G4float  GREfactor = 5.0f;
const G4float  GRWfactor = 0.3f;
const G4double MinGammaEnergy = 10*CLHEP::keV;
const G4double MaxDeltaEnergy = CLHEP::MeV;
const G4double LevelDensity   = 0.125/CLHEP::MeV;
const G4double NormC = 2.5*CLHEP::millibarn/(CLHEP::pi2*CLHEP::hbarc*CLHEP::hbarc);
const G4double Tolerance = 0.1*CLHEP::keV;

G4PhotonEvaporation::G4PhotonEvaporation(G4GammaTransition* p)
  : fLevelManager(nullptr), fTransition(p), fVerbose(0), fPoints(0), 
    vShellNumber(-1), fIndex(0), fTimeLimit(DBL_MAX), fMaxLifeTime(DBL_MAX), 
    fICM(false), fRDM(false), fSampleTime(true)
{
  fNuclearLevelData = G4NuclearLevelData::GetInstance(); 
  if(!fTransition) { 
    char* en = getenv("G4UseNuclearPolarization"); 
    if(en) { fTransition = new G4PolarizedGammaTransition(); }
    else   { fTransition = new G4GammaTransition(); }
  }

  char* env = getenv("G4AddTimeLimitToPhotonEvaporation"); 
  if(env) { fTimeLimit = 1.e-16*CLHEP::second; }

  theA = theZ = fCode = 0;
  fLevelEnergyMax = fStep = fExcEnergy = fFermiEnergy = fProbability = 0.0;

  for(G4int i=0; i<MAXDEPOINT; ++i) { fCummProbability[i] = 0.0; }
  if(0.0f == GREnergy[1]) {
    G4Pow* g4pow = G4Pow::GetInstance();
    for (G4int A=1; A<MAXGRDATA; ++A) {
      GREnergy[A] = (G4float)(40.3*CLHEP::MeV/g4pow->powZ(A,0.2));
      GRWidth[A] = GRWfactor*GREnergy[A];
    }
  } 
}

G4PhotonEvaporation::~G4PhotonEvaporation()
{ 
  delete fTransition;
}

G4Fragment* 
G4PhotonEvaporation::EmittedFragment(G4Fragment* nucleus)
{
  if(fRDM) { fSampleTime = false; }
  else     { fSampleTime = true; }

  G4Fragment* gamma = GenerateGamma(nucleus);
  if(fVerbose > 0) {
    G4cout << "G4PhotonEvaporation::EmittedFragment: RDM= " << fRDM << G4endl;   
    if(gamma) { G4cout << *gamma << G4endl; }
    G4cout << "   Residual: " << *nucleus << G4endl;
  }
  return gamma; 
}

G4FragmentVector* 
G4PhotonEvaporation::BreakUpFragment(G4Fragment* nucleus)
{
  //G4cout << "G4PhotonEvaporation::BreakUpFragment" << G4endl;
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, nucleus);
  return products;
}

G4FragmentVector* 
G4PhotonEvaporation::BreakUp(const G4Fragment& nucleus)
{
  //G4cout << "G4PhotonEvaporation::BreakUp" << G4endl;
  G4Fragment* aNucleus = new G4Fragment(nucleus);
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, aNucleus);
  products->push_back(aNucleus);
  return products;
}

G4FragmentVector* 
G4PhotonEvaporation::BreakItUp(const G4Fragment& nucleus)
{
  //G4cout << "G4PhotonEvaporation::BreakItUp" << G4endl;
  G4Fragment* aNucleus = new G4Fragment(nucleus);
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, aNucleus);
  products->push_back(aNucleus);
  return products;
}

G4bool G4PhotonEvaporation::BreakUpChain(G4FragmentVector* products,
					 G4Fragment* nucleus)
{
  if(fVerbose > 0) {
    G4cout << "G4PhotonEvaporation::BreakUpChain RDM= " << fRDM << " "
	   << *nucleus << G4endl;
  }
  G4Fragment* gamma = 0;
  if(fRDM) { fSampleTime = false; }
  else     { fSampleTime = true; }

  do {
    gamma = GenerateGamma(nucleus);
    if(gamma) { 
      products->push_back(gamma);
      if(fVerbose > 0) {
	G4cout << "G4PhotonEvaporation::BreakUpChain: "   
	       << *gamma << G4endl;
	G4cout << "   Residual: " << *nucleus << G4endl;
      }
      // for next decays in the chain always sample time
      fSampleTime = true;
    } 
    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  } while(gamma);
  return false;
}

G4double 
G4PhotonEvaporation::GetEmissionProbability(G4Fragment* nucleus) 
{
  fProbability = 0.0;
  fExcEnergy = nucleus->GetExcitationEnergy();
  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  fCode = 1000*Z + A;
  if(fVerbose > 1) {
    G4cout << "G4PhotonEvaporation::GetEmissionProbability: Z=" 
	   << Z << " A=" << A << " Eexc(MeV)= " << fExcEnergy << G4endl; 
  }
  // ignore gamma de-excitation for exotic fragments 
  // and for very low excitations
  if(0 >= Z || 1 >= A || Z == A || Tolerance >= fExcEnergy)
    { return fProbability; }

  // ignore gamma de-excitation for highly excited levels
  if(A >= MAXGRDATA) { A =  MAXGRDATA-1; }
  //G4cout << " GREnergy= " << GREnergy[A] << " GRWidth= " << GRWidth[A] << G4endl; 
  if(fExcEnergy >= (G4double)(GREfactor*GRWidth[A] + GREnergy[A])) { 
    return fProbability; 
  }

  // probability computed assuming continium transitions
  // VI: continium transition are limited only to final states
  //     below Fermi energy (this approach needs further evaluation)
  fFermiEnergy = std::max(0.0, nucleus->ComputeGroundStateMass(Z, A-1) 
    + CLHEP::neutron_mass_c2 - nucleus->GetGroundStateMass());

  // max energy level for continues transition
  G4double emax = std::min(fFermiEnergy, fExcEnergy);
  const G4double eexcfac = 0.99;
  if(0.0 == emax || fExcEnergy*eexcfac <= emax) { emax = fExcEnergy*eexcfac; }

  fStep = emax;
  fPoints = std::min((G4int)(fStep/MaxDeltaEnergy) + 2, MAXDEPOINT);
  fStep /= ((G4double)(fPoints - 1));
  if(fVerbose > 1) {
    G4cout << "Emax= " << emax << " Npoints= " << fPoints 
	   <<" Efermi= " << fFermiEnergy << "  Eex= " << fExcEnergy << G4endl;
  }
  // integrate probabilities
  G4double eres = (G4double)GREnergy[A];
  G4double wres = (G4double)GRWidth[A];
  G4double eres2= eres*eres;
  G4double wres2= wres*wres;
  G4double xsqr = std::sqrt(A*LevelDensity*fExcEnergy);

  G4double egam    = fExcEnergy - emax;
  G4double gammaE2 = egam*egam;
  G4double gammaR2 = gammaE2*wres*wres;
  G4double egdp2   = gammaE2 - eres*eres;

  G4double p0 = G4Exp(2*(std::sqrt(A*LevelDensity*emax) - xsqr))
    *gammaR2*gammaE2/(egdp2*egdp2 + gammaR2);
  G4double p1(0.0);

  for(G4int i=1; i<fPoints; ++i) {
    egam += fStep;
    gammaE2 = egam*egam;
    gammaR2 = gammaE2*wres2;
    egdp2   = gammaE2 - eres2;
    //G4cout << "Egamma= " << egam << "  Eex= " << fExcEnergy
    //<< "  p0= " << p0 << " p1= " << p1 << G4endl;
    p1 = G4Exp(2*(std::sqrt(A*LevelDensity*std::abs(fExcEnergy - egam)) - xsqr))
      *gammaR2*gammaE2/(egdp2*egdp2 + gammaR2);
    fProbability += (p1 + p0);
    fCummProbability[i] = fProbability;
    p0 = p1;
  }
  fProbability *= 0.5*fStep*NormC*A;
  if(fVerbose > 1) { G4cout << "prob= " << fProbability << G4endl; }
  return fProbability;
}

G4double 
G4PhotonEvaporation::GetFinalLevelEnergy(G4int Z, G4int A, G4double energy)
{
  G4double E = energy;
  InitialiseLevelManager(Z, A);
  if(fLevelManager) { 
    E = (G4double)fLevelManager->NearestLevelEnergy(energy, fIndex); 
    if(E > fLevelEnergyMax + Tolerance) { E = energy; }
  }
  return E;
}

G4double G4PhotonEvaporation::GetUpperLevelEnergy(G4int Z, G4int A)
{
  InitialiseLevelManager(Z, A);
  return fLevelEnergyMax;
}

G4Fragment* 
G4PhotonEvaporation::GenerateGamma(G4Fragment* nucleus)
{
  G4Fragment* result = 0;
  G4double eexc = nucleus->GetExcitationEnergy();
  if(eexc < MinGammaEnergy) { return result; }

  InitialiseLevelManager(nucleus->GetZ_asInt(), nucleus->GetA_asInt());

  G4double time = nucleus->GetCreationTime();

  G4double efinal = 0.0;
  vShellNumber    = -1;
  size_t shell    = 0;
  G4int  deltaS   = 1;
  G4bool isGamma  = true;
  G4bool isLongLived  = false;
  G4bool isX = false;
  G4bool icm = fICM;

  if(fVerbose > 1) {
    G4cout << "GenerateGamma: Exc= " << eexc << " Emax= " 
	   << fLevelEnergyMax << G4endl;
  }
  // continues part
  if(!fLevelManager || eexc > fLevelEnergyMax + Tolerance) {
    //G4cout << "Continues fPoints= " << fPoints << " " << fLevelManager << G4endl;

    // we compare current excitation versus value used for probability 
    // computation and also Z and A used for probability computation 
    if(fCode != 1000*theZ + theA && eexc != fExcEnergy) { 
      GetEmissionProbability(nucleus); 
    }
    if(fProbability == 0.0) { return result; }
    G4double y = fProbability*G4UniformRand();
    for(G4int i=1; i<fPoints; ++i) {
      //G4cout << "y= " << y << " cummProb= " << fCummProbability[i] << G4endl;
      if(y <= fCummProbability[i]) {
	efinal = fStep*((i - 1) + (y - fCummProbability[i-1])
			/(fCummProbability[i] - fCummProbability[i-1]));
	break;
      }
    }
    // final discrete level
    if(fLevelManager && efinal <=  fLevelEnergyMax + Tolerance) {
      //G4cout << "Efinal= " << efinal << "  idx= " << fIndex << G4endl;
      fIndex = fLevelManager->NearestLevelIndex(efinal, fIndex);
      efinal = (G4double)fLevelManager->LevelEnergy(fIndex);
    }
    //discrete part
  } else { 
    fIndex = fLevelManager->NearestLevelIndex(eexc, fIndex);
    if(fVerbose > 1) {
      G4cout << "Discrete emission from level Index= " << fIndex 
	     << " Elevel= " << fLevelManager->LevelEnergy(fIndex)
	     << "  ICM= " << fICM << G4endl;
    }
    if(0 == fIndex) { return result; }
    G4double ltime = 0.0;
    if(fICM) { ltime = (G4double)fLevelManager->LifeTime(fIndex); }
    else     { ltime = (G4double)fLevelManager->LifeTimeGamma(fIndex); }

    if(ltime >= fMaxLifeTime) { return result; }
    if(ltime > fTimeLimit) { 
      icm = true; 
      isLongLived = true;
    }
    const G4NucLevel* level = fLevelManager->GetLevel(fIndex);
    size_t ntrans = level->NumberOfTransitions();
    size_t idx = 0;
    isX = level->IsXLevel();
    //G4cout << "Ntrans= " << ntrans << " idx= " << idx << " isX= " << isX 
    //	   << " icm= " << icm << G4endl;
    if(!isX) {
      if(1 < ntrans) {
	G4double rndm = G4UniformRand();
	if(icm) { idx = level->SampleGammaETransition(rndm); }
	else    { idx = level->SampleGammaTransition(rndm); }
	//G4cout << "Sampled idx= " << idx << "  rndm= " << rndm << G4endl;
      }
      if(icm) {
	G4double rndm = G4UniformRand();
	G4double prob = level->GammaProbability(idx);
	if(rndm > prob) {
	  rndm = (rndm - prob)/(1.0 - prob);
	  shell = level->SampleShell(idx, rndm);
	  vShellNumber = shell;
          isGamma = false;
	}
      }
    }
    efinal = (G4double)level->FinalExcitationEnergy(idx);
    fIndex = idx;

    if(fSampleTime && ltime > 0.0) { 
      time -= ltime*G4Log(G4UniformRand()); 
    }
  }
  // no gamma emission for X-level
  if(isX) {
    G4ThreeVector v = nucleus->GetMomentum().vect().unit();
    G4double mass = nucleus->GetGroundStateMass() + efinal;
    G4double e = std::max(mass,nucleus->GetMomentum().e()); 
    G4double mom = std::sqrt((e - mass)*(e + mass)); 
    v *= mom;
    nucleus->SetMomentum(G4LorentzVector(v.x(),v.y(),v.z(),efinal+mass));
    // normal case
  } else { 
    result = fTransition->SampleTransition(nucleus, efinal,
					   deltaS, shell,
					   isGamma, isLongLived);
    if(result) { result->SetCreationTime(time); }
  }
  nucleus->SetCreationTime(time);
  
  if(fVerbose > 1) { 
    G4cout << "Final level E= " << efinal << " time= " << time 
	   << " idx= " << fIndex << " isX " << isX 
	   << " isGamma: " << isGamma << " isLongLived: " << isLongLived
	   << " deltaS= " << deltaS << " shell= " << shell << G4endl;
  }
  return result;
}

void G4PhotonEvaporation::SetMaxHalfLife(G4double val)
{
  static const G4double tfact = G4Pow::GetInstance()->logZ(2);
  fMaxLifeTime = val/tfact;
}

void G4PhotonEvaporation::SetGammaTransition(G4GammaTransition* p)
{
  if(p != fTransition) {
    delete fTransition;
    fTransition = p;
  }
}

void G4PhotonEvaporation::SetICM(G4bool val)
{
  fICM = val;
}

void G4PhotonEvaporation::RDMForced(G4bool val)
{
  fRDM = val;
}
  
