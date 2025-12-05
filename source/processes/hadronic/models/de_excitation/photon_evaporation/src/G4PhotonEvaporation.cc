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
// -------------------------------------------------------------------
//

#include "G4PhotonEvaporation.hh"

#include "G4NuclearPolarizationStore.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"
#include "G4FragmentVector.hh"
#include "G4GammaTransition.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsModelCatalog.hh"

G4float G4PhotonEvaporation::GREnergy[] = {0.0f};
G4float G4PhotonEvaporation::GRWidth[] = {0.0f};

namespace
{
  constexpr G4double timeLimit = 10*CLHEP::ns;
  constexpr G4double eLimit = 200*CLHEP::keV;
}

G4PhotonEvaporation::G4PhotonEvaporation(G4GammaTransition* p)
  : fTransition(p), fPolarization(nullptr), fVerbose(1)
{
  if (fVerbose > 1) {
    G4cout << "### New G4PhotonEvaporation() " << this << G4endl;
  }
  fNuclearLevelData = G4NuclearLevelData::GetInstance(); 
  fTolerance = 20*CLHEP::eV;
  fCummProbability[0] = 0.0;
  if(nullptr == fTransition) { fTransition = new G4GammaTransition(); }

  fSecID = G4PhysicsModelCatalog::GetModelID("model_G4PhotonEvaporation");

  if(0.0f == GREnergy[2]) { InitialiseGRData(); }
}

G4PhotonEvaporation::~G4PhotonEvaporation()
{ 
  delete fTransition;
}

void G4PhotonEvaporation::Initialise()
{
  if (isInitialised) { return; }
  isInitialised = true;

  G4DeexPrecoParameters* param = fNuclearLevelData->GetParameters();
  fTolerance = param->GetMinExcitation();
  fMaxLifeTime = param->GetMaxLifeTime();
  fLocalTimeLimit = fRDM ? fMaxLifeTime : std::max(fMaxLifeTime, timeLimit);
  fCorrelatedGamma = param->CorrelatedGamma();
  fICM = param->GetInternalConversionFlag();
  fVerbose = param->GetVerbose();

  fTransition->SetPolarizationFlag(fCorrelatedGamma);
  fTransition->SetTwoJMAX(param->GetTwoJMAX());
  fTransition->SetVerbose(fVerbose);
  if (fVerbose > 1) {
    G4cout << "### G4PhotonEvaporation is initialized " << this << G4endl;   
  }
}

void G4PhotonEvaporation::InitialiseGRData()
{
  if(0.0f == GREnergy[2]) {
    G4Pow* g4calc = G4Pow::GetInstance();
    const G4float GRWfactor = 0.3f;
    for (G4int A=1; A<MAXGRDATA; ++A) {
      GREnergy[A] = (G4float)(40.3*CLHEP::MeV/g4calc->powZ(A,0.2));
      GRWidth[A] = GRWfactor*GREnergy[A];
    }
  }
}

G4Fragment* 
G4PhotonEvaporation::EmittedFragment(G4Fragment* nucleus)
{
  if(!isInitialised) { Initialise(); }
  fSampleTime = !fRDM;
  if (fRDM) {
    nucleus->SetNumberOfElectrons(nucleus->GetZ_asInt());
  }

  // potentially external code may set initial polarization
  // but only for radioactive decay nuclear polarization is considered
  G4NuclearPolarizationStore* fNucPStore = nullptr;
  if(fCorrelatedGamma && fRDM) {
    fNucPStore = G4NuclearPolarizationStore::GetInstance();
    auto nucp = nucleus->GetNuclearPolarization();
    if(nullptr != nucp) { 
      fNucPStore->RemoveMe(nucp);
    } 
    fPolarization = fNucPStore->FindOrBuild(nucleus->GetZ_asInt(),
                                            nucleus->GetA_asInt(),
                                            nucleus->GetExcitationEnergy());
    nucleus->SetNuclearPolarization(fPolarization);
  }
  if(fVerbose > 2) { 
    G4cout << "G4PhotonEvaporation::EmittedFragment: " 
           << *nucleus << G4endl;
    if (nullptr != fPolarization) { G4cout << "NucPolar: " << fPolarization << G4endl; }
    G4cout << " CorrGamma: " << fCorrelatedGamma << " RDM: " << fRDM
           << " fPolarization: " << fPolarization << G4endl;
  }
  G4Fragment* gamma = GenerateGamma(nucleus);
  
  if(gamma != nullptr) { gamma->SetCreatorModelID(fSecID); }
  
  // remove G4NuclearPolarizaton when reach ground state
  if (nullptr != fNucPStore && nullptr != fPolarization && 0 == fIndex) {
    if(fVerbose > 3) { 
      G4cout << "G4PhotonEvaporation::EmittedFragment: remove " 
             << fPolarization << G4endl;
    }
    fNucPStore->RemoveMe(fPolarization);
    fPolarization = nullptr;
    nucleus->SetNuclearPolarization(fPolarization);
  }

  if(fVerbose > 2) {
    G4cout << "G4PhotonEvaporation::EmittedFragment: RDM= " 
           << fRDM << " done:" << G4endl; 
    if(gamma) { G4cout << *gamma << G4endl; }
    G4cout << "   Residual: " << *nucleus << G4endl;
  }
  return gamma; 
}

G4FragmentVector* 
G4PhotonEvaporation::BreakItUp(const G4Fragment& nucleus)
{
  if(fVerbose > 1) {
    G4cout << "G4PhotonEvaporation::BreakItUp" << G4endl;
  }
  G4Fragment* aNucleus = new G4Fragment(nucleus);
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, aNucleus);
  aNucleus->SetCreatorModelID(fSecID);
  products->push_back(aNucleus);
  return products;
}

G4bool G4PhotonEvaporation::BreakUpChain(G4FragmentVector* products,
                                         G4Fragment* nucleus)
{
  if(!isInitialised) { Initialise(); }
  if(fVerbose > 1) {
    G4cout << "G4PhotonEvaporation::BreakUpChain RDM= " << fRDM << " "
           << *nucleus << G4endl;
  }
  G4Fragment* gamma = nullptr;
  fSampleTime = !fRDM;

  // start decay chain from unpolarized state
  if(fCorrelatedGamma) {
    fPolarization = new G4NuclearPolarization(nucleus->GetZ_asInt(),
                                              nucleus->GetA_asInt(),
                                              nucleus->GetExcitationEnergy());
    nucleus->SetNuclearPolarization(fPolarization);
  }

  do {
    gamma = GenerateGamma(nucleus);
    if (nullptr != gamma) {
      gamma->SetCreatorModelID(fSecID);
      products->push_back(gamma);
    }
    // for next decays in the chain always sample time
    fSampleTime = true;
    if (fVerbose > 2) {
      G4cout << "G4PhotonEvaporation::BreakUpChain: next decay" << G4endl;
      if (nullptr != gamma) { G4cout << "   " << *gamma << G4endl; }
      else { G4cout << "   not possible" << G4endl; }
      G4cout << "   Residual: " << *nucleus << G4endl;
    }
    // Loop checking, 22-Dec-2024, Vladimir Ivanchenko
  } while (!(nucleus->IsLongLived() || nucleus->GetExcitationEnergy() <= fTolerance));

  // clear nuclear polarization end of chain
  if(nullptr != fPolarization) {
    delete fPolarization;
    fPolarization = nullptr;
    nucleus->SetNuclearPolarization(fPolarization);
  }  
  return false;
}

G4double 
G4PhotonEvaporation::GetEmissionProbability(G4Fragment* nucleus) 
{
  if(!isInitialised) { Initialise(); }
  fProbability = 0.0;
  fExcEnergy = nucleus->GetExcitationEnergy();
  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  if(fVerbose > 2) {
    G4cout << "G4PhotonEvaporation::GetEmissionProbability: Z=" 
           << Z << " A=" << A << " Eexc(MeV)= " << fExcEnergy << G4endl; 
  }
  // ignore gamma de-excitation for exotic fragments 
  // and for very low excitations
  if(0 >= Z || 1 >= A || Z == A || fTolerance >= fExcEnergy)
    { return fProbability; }

  // ignore gamma de-excitation for highly excited levels
  if(A >= MAXGRDATA) { A =  MAXGRDATA-1; }

  static const G4double GREfactor = 5.0;
  G4double edelta = GREfactor*(G4double)GRWidth[A] + (G4double)GREnergy[A];
  if (fVerbose > 2)
    G4cout << "   GREnergy=" << GREnergy[A] << " GRWidth="<<GRWidth[A]
	   << " Edelta=" << edelta <<G4endl; 
  if (fExcEnergy >= edelta) { 
    return fProbability; 
  }

  // probability computed assuming continium transitions in the frame of the nucleus
  fStep = fExcEnergy;
  const G4double MaxDeltaEnergy = CLHEP::MeV;
  fPoints = std::min((G4int)(fStep/MaxDeltaEnergy) + 2, MAXDEPOINT);
  fStep /= ((G4double)(fPoints - 1));

  if(fVerbose > 2) {
    G4cout << "  Npoints= " << fPoints 
           << "  Eex=" << fExcEnergy << " Estep=" << fStep << G4endl;
  }

  // integrate probabilities
  G4double eres = (G4double)GREnergy[A];
  G4double wres = (G4double)GRWidth[A];
  G4double eres2= eres*eres;
  G4double wres2= wres*wres;

  // initial state
  G4double levelDensity = fNuclearLevelData->GetLevelDensity(Z,A,fExcEnergy);
  G4double xdrt = G4Exp(2*std::sqrt(levelDensity*fExcEnergy));

  // the loop over excitation energy of the residual nucleus
  // from 0 to fExcEnergy
  // gamma energy is defined via non-relativistic formula
  G4double egam    = fExcEnergy;
  G4double gammaE2 = egam*egam;
  G4double gammaR2 = gammaE2*wres2;
  G4double egdp2   = gammaE2 - eres2;

  G4double p0 = egam*gammaR2*gammaE2/(egdp2*egdp2 + gammaR2);
  G4double p1, e;

  for(G4int i=1; i<fPoints; ++i) {
    egam -= fStep;
    if (i + 1 == fPoints) {
      p1 = 0.0;
    } else {
      gammaE2 = egam*egam;
      gammaR2 = gammaE2*wres2;
      egdp2   = gammaE2 - eres2;
      e = fExcEnergy - egam;
      levelDensity = fNuclearLevelData->GetLevelDensity(Z, A, e);
      p1 = egam*G4Exp(2.0*(std::sqrt(levelDensity*e)))*gammaR2*gammaE2/(egdp2*egdp2 + gammaR2);
    }
    fProbability += (p1 + p0);
    fCummProbability[i] = fProbability;
    if(fVerbose > 3) {
      G4cout << "Egamma= " << egam << "  Eex= " << fExcEnergy
             << "  p0= " << p0 << " p1= " << p1 << " sum= " 
             << fCummProbability[i] <<G4endl;
    }
    p0 = p1;
  }

  static const G4double NormC = 1.25*CLHEP::millibarn
    /(CLHEP::pi2*CLHEP::hbarc*CLHEP::hbarc);
  fProbability *= fStep*NormC*A/xdrt;
  if(fVerbose > 1) { G4cout << "prob= " << fProbability << G4endl; }
  return fProbability;
}

G4double
G4PhotonEvaporation::ComputeInverseXSection(G4Fragment*, G4double)
{
  return 0.0;
}

G4double 
G4PhotonEvaporation::ComputeProbability(G4Fragment* theNucleus, G4double)
{
  return GetEmissionProbability(theNucleus);
}

G4double
G4PhotonEvaporation::GetFinalLevelEnergy(G4int Z, G4int A, G4double energy)
{
  G4double E = energy;
  InitialiseLevelManager(Z, A);
  if (nullptr != fLevelManager) { 
    E = fLevelManager->NearestLevelEnergy(energy, fIndex); 
    if(E > fLevelEnergyMax + fTolerance) { E = energy; }
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
  if(!isInitialised) { Initialise(); }
  G4Fragment* result = nullptr;

  // initial level
  G4double eexc = nucleus->GetExcitationEnergy();
  InitialiseLevelManager(nucleus->GetZ_asInt(), nucleus->GetA_asInt());
  // long life time flag - "true" for a fragment, which will be tracked
  G4bool isLL = false;
  // lifetime of the fragment
  G4double ltime = 0.0;
  fExcEnergy = eexc;
  // index is unknown - default is the ground state
  fIndex = 0;

  G4double time = nucleus->GetCreationTime();
  G4double elevel = eexc;
  G4double efinal = 0.0;
  G4double ratio  = 0.0;
  vShellNumber    = -1;
  G4int  JP1      = 0;
  G4int  JP2      = 0;
  G4int  multiP   = 0;
  G4bool isGamma  = true;
  G4bool isDiscrete = false;
  G4bool finalDiscrete = false;

  const G4NucLevel* level = nullptr;
  std::size_t ntrans = 0;

  if(fVerbose > 2) {
    G4cout << "## GenerateGamma: Z=" << theZ << " A=" << theA << " Eex= " << eexc
           << " Eexmax= " << fLevelEnergyMax << G4endl;
  }
  // initial discrete state is ground level
  if (eexc <= fTolerance) {
    isDiscrete = true;

    // initial state may be a discrete level
  } else if (nullptr != fLevelManager && eexc <= fLevelEnergyMax + fTolerance) {
    fIndex = fLevelManager->NearestLevelIndex(eexc);
    elevel = fLevelManager->LevelEnergy(fIndex);
    isDiscrete = (std::abs(elevel - eexc) < fTolerance);
    if(fVerbose > 2) {
      G4cout << "          Level index=" << fIndex 
             << " lTime=" << fLevelManager->LifeTime(fIndex)
	     << " Elevel=" << elevel
	     << " isDiscrete:" << isDiscrete << G4endl;
    }
    if(isDiscrete && 0 < fIndex) {
      // for discrete transition  
      level = fLevelManager->GetLevel(fIndex);
      if(nullptr != level) { 
        ntrans = level->NumberOfTransitions();
	G4int idxfl = fLevelManager->FloatingLevel(fIndex); 
	// for floating level check levels with the same energy
        if (idxfl > 0) {
	  auto newlevel = fLevelManager->GetLevel(fIndex-1);
	  G4double newenergy = fLevelManager->LevelEnergy(fIndex-1);
	  if (nullptr != newlevel && std::abs(elevel - newenergy) < fTolerance) {
	    std::size_t newntrans = newlevel->NumberOfTransitions();
	    if (newntrans > 0) {
	      --fIndex;
	      level = newlevel;
	      elevel = newenergy;
	      ntrans = newntrans;
	    }
	  }
	}
        JP1 = std::abs(fLevelManager->TwoSpinParity(fIndex)); 
        if(fVerbose > 2) {
          G4cout << "          ntrans= " << ntrans << " JP= " << JP1
                 << " RDM: " << fRDM << G4endl;
        }
      }
      // if a level has no defined transitions
      if (0 == ntrans) {
	isDiscrete = false;
      }
      // transition from continues spectrum to the ground state
    } else if (0 == fIndex) {
      isDiscrete = true;
    }
  }

  if(fVerbose > 2) {
    G4long prec = G4cout.precision(4);
    G4cout << "   Z=" << nucleus->GetZ_asInt()
           << " A=" << nucleus->GetA_asInt() 
           << " Exc=" << eexc << " Emax=" 
           << fLevelEnergyMax << " idx=" << fIndex
           << " fPoints= " << fPoints
           << " Ntr=" << ntrans << " discrete:" << isDiscrete
           << G4endl;
    G4cout.precision(prec);
  }

  if(!isDiscrete) {
    // primary fragment is in continium 
    GetEmissionProbability(nucleus); 

    if(fProbability == 0.0) { 
      fPoints = 1; 
      efinal = 0.0; 
    } else {
      G4double y = fCummProbability[fPoints-1]*G4UniformRand();
      for(G4int i=1; i<fPoints; ++i) {
        if(fVerbose > 3) {
          G4cout << "y= " << y << " cummProb= " << fCummProbability[i] 
                 << " fPoints= " << fPoints << " fStep= " << fStep << G4endl;
        }
        if(y <= fCummProbability[i]) {
          efinal = fStep*((i - 1) + (y - fCummProbability[i-1])
                          /(fCummProbability[i] - fCummProbability[i-1]));
          break;
        }
      }
    }
    // final discrete level or continues exitation energy
    if(fVerbose > 2) {
      G4cout << "Continues proposes Efinal=" << efinal
	     << " Initial Idx=" << fIndex << G4endl;
    }
    
    if(nullptr != fLevelManager) {
      // final discrete level
      if (efinal < fLevelEnergyMax + fTolerance) {
        fIndex = fLevelManager->NearestLevelIndex(efinal, fIndex);
        G4double el = fLevelManager->LevelEnergy(fIndex);
	// protection - take level below 
	if (el >= eexc + fTolerance && 0 < fIndex) {
	  --fIndex;
	  el = fLevelManager->LevelEnergy(fIndex);
	}
	// further decays will be discrete
        ltime = fLevelManager->LifeTime(fIndex);
	if (fIndex <= 1 || std::abs(efinal - el) <= eLimit || ltime >= fLocalTimeLimit) {
	  efinal = el;
	  finalDiscrete = true;
	} else {
	  fIndex = 0;
	}
      }
    }
    if (fVerbose > 2) {
      G4cout << "Continues emission efinal(MeV)= " << efinal
	     << " idxFinal=" << fIndex << " isdiscrete:" << isDiscrete << G4endl; 
    }

    // initial continues and final ground state
  } else if (0 == fIndex) {
    efinal = 0.0;
    isDiscrete = false;
    if (nullptr != fLevelManager) { finalDiscrete = true; }

    // discrete part for excited nucleus
  } else {
 
    if (fVerbose > 2) {
      G4cout << "Discrete emission from level Index=" << fIndex 
             << " Elevel=" << fLevelManager->LevelEnergy(fIndex)
             << " Ltime=" << fLevelManager->LifeTime(fIndex)
             << " LtimeMax=" << fLocalTimeLimit
             << "  RDM=" << fRDM << " ICM=" << fICM << G4endl;
    }

    // stable fragment has life time DBL_MAX
    ltime = fLevelManager->LifeTime(fIndex);

    // stable isomer - no sampling of transition 
    if (ltime == DBL_MAX) {
      nucleus->SetFloatingLevelNumber(0);
      nucleus->SetLongLived(true);
      return result;
    }

    // sampling index of a final level
    std::size_t idx = 0;
    if(1 < ntrans) {
      idx = level->SampleGammaTransition(G4UniformRand());
    }
    if(fVerbose > 2) {
      G4cout << "Ntrans= " << ntrans << " idx= " << idx
	     << " ICM= " << fICM << "  abs(JP1)= " << JP1 << G4endl;
    }

    // sampling IC or gamma transition
    G4double prob = (G4double)level->GammaProbability(idx);

    // prob = 0 means that there is only internal conversion
    if (prob < 1.0) {
      G4double rndm = G4UniformRand();
      if (rndm > prob) {
	isGamma = false;
	if (fICM) {
	  rndm = (rndm - prob)/(1.0 - prob);
	  vShellNumber = level->SampleShell(idx, rndm);
	}
      }
    }
    // it is a discrete transition with possible gamma correlation
    ratio  = level->MultipolarityRatio(idx);
    multiP = level->TransitionType(idx);
    fIndex = level->FinalExcitationIndex(idx);
    finalDiscrete = true;

    // final level parameters
    efinal = fLevelManager->LevelEnergy(fIndex);
    // time is sampled if decay not prompt and this class called not 
    // from radioactive decay and isomer production is enabled 
    if(fSampleTime && ltime > 0.0) {
      time -= ltime*G4Log(G4UniformRand()); 
    }
  }
  ltime = 0.0;
  if (finalDiscrete) {
    ltime = fLevelManager->LifeTime(fIndex);
    JP2 = fLevelManager->TwoSpinParity(fIndex);
  }
  // sample continues or discrete transition if transition
  // is above distance between floating level
  if (std::abs(efinal - eexc) > fTolerance) {
    result = fTransition->SampleTransition(nucleus, efinal, ratio, JP1,
					   std::abs(JP2), multiP, vShellNumber, 
					   isDiscrete, isGamma);
    if (nullptr != result) { result->SetCreationTime(time); }
  }
  // update parameters of the fragment
  nucleus->SetCreationTime(time);
  nucleus->SetSpin(0.5*JP2);
  if (nullptr != fPolarization) { fPolarization->SetExcitationEnergy(efinal); }
    
  if (finalDiscrete) {
    G4int idxfl = fLevelManager->FloatingLevel(fIndex);
    nucleus->SetFloatingLevelNumber(idxfl);

    if (ltime > fLocalTimeLimit) { isLL = true; }
  }
  nucleus->SetLongLived(isLL);
      
  if (fVerbose > 2) {
    G4String ss = "## ";
    if (isLL && efinal > 0.0 && efinal < MeV) { ss += "=I="; } 
    if (isLL && efinal >= MeV) { ss += "=J="; } 
    if (efinal >= 6*MeV) { ss += "=K="; } 
    G4cout << "   " << ss << " Efinal=" << efinal 
	   << " Efrag=" << nucleus->GetExcitationEnergy()
           << " lt=" << ltime 
           << " idxFin=" << fIndex << " isDiscrete:" << isDiscrete
           << " isGamma:" << isGamma << " isStable:" << isLL
	   << " multiP=" << multiP << " shell=" << vShellNumber 
           << " abs(JP1)= " << JP1 << " abs(JP2)= " << JP2 << G4endl;
  }
  return result;
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
  
