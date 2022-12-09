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
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include "G4PhysicsModelCatalog.hh"
#include "G4AutoLock.hh"

namespace
{
  G4Mutex photEvaporationMutex = G4MUTEX_INITIALIZER;
}

G4float G4PhotonEvaporation::GREnergy[] = {0.0f};
G4float G4PhotonEvaporation::GRWidth[] = {0.0f};

G4PhotonEvaporation::G4PhotonEvaporation(G4GammaTransition* p)
  : fLevelManager(nullptr), fTransition(p), fPolarization(nullptr),
    fVerbose(1), fPoints(0), vShellNumber(-1), fIndex(0), fSecID(-1),
    fMaxLifeTime(DBL_MAX), 
    fICM(true), fRDM(false), fSampleTime(true), 
    fCorrelatedGamma(false), isInitialised(false)
{
  //G4cout << "### New G4PhotonEvaporation() " << this << G4endl;   
  fNuclearLevelData = G4NuclearLevelData::GetInstance(); 
  fTolerance = 20*CLHEP::eV;

  if(!fTransition) { fTransition = new G4GammaTransition(); }

  theA = theZ = fCode = 0;
  fSecID = G4PhysicsModelCatalog::GetModelID("model_G4PhotonEvaporation");
  fLevelEnergyMax = fStep = fExcEnergy = fProbability = 0.0;

  for(G4int i=0; i<MAXDEPOINT; ++i) { fCummProbability[i] = 0.0; }
  if(0.0f == GREnergy[1]) { InitialiseGRData(); }
}

G4PhotonEvaporation::~G4PhotonEvaporation()
{ 
  delete fTransition;
}

void G4PhotonEvaporation::Initialise()
{
  if(isInitialised) { return; }
  isInitialised = true;

  G4DeexPrecoParameters* param = fNuclearLevelData->GetParameters();
  fTolerance = param->GetMinExcitation();
  fMaxLifeTime = param->GetMaxLifeTime();
  fCorrelatedGamma = param->CorrelatedGamma();
  fICM = param->GetInternalConversionFlag();
  fVerbose = param->GetVerbose();

  fTransition->SetPolarizationFlag(fCorrelatedGamma);
  fTransition->SetTwoJMAX(param->GetTwoJMAX());
  fTransition->SetVerbose(fVerbose);
  if(fVerbose > 1) {
    G4cout << "### G4PhotonEvaporation is initialized " << this << G4endl;   
  }
}

void G4PhotonEvaporation::InitialiseGRData()
{
  if(0.0f == GREnergy[1]) {
    G4AutoLock l(&photEvaporationMutex);
    if(0.0f == GREnergy[1]) {
      G4Pow* g4calc = G4Pow::GetInstance();
      const G4float GRWfactor = 0.3f;
      for (G4int A=1; A<MAXGRDATA; ++A) {
	GREnergy[A] = (G4float)(40.3*CLHEP::MeV/g4calc->powZ(A,0.2));
	GRWidth[A] = GRWfactor*GREnergy[A];
      }
    }
    l.unlock();    
  }
}

G4Fragment* 
G4PhotonEvaporation::EmittedFragment(G4Fragment* nucleus)
{
  if(!isInitialised) { Initialise(); }
  fSampleTime = !fRDM;

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
    if(fPolarization) { G4cout << "NucPolar: " << fPolarization << G4endl; }
    G4cout << " CorrGamma: " << fCorrelatedGamma << " RDM: " << fRDM
           << " fPolarization: " << fPolarization << G4endl;
  }
  G4Fragment* gamma = GenerateGamma(nucleus);
  
  if(gamma != nullptr) { gamma->SetCreatorModelID(fSecID); }
  
  // remove G4NuclearPolarizaton when reach ground state
  if(fNucPStore && fPolarization && 0 == fIndex) {
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
    if(gamma) {
      gamma->SetCreatorModelID(fSecID);
      products->push_back(gamma);
      if(fVerbose > 2) {
        G4cout << "G4PhotonEvaporation::BreakUpChain: "   
               << *gamma << G4endl;
        G4cout << "   Residual: " << *nucleus << G4endl;
      }
      // for next decays in the chain always sample time
      fSampleTime = true;
    } 
    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  } while(gamma);

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
  fCode   = 1000*Z + A; 
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
  //G4cout<<" GREnergy= "<< GREnergy[A]<<" GRWidth= "<<GRWidth[A]<<G4endl; 

  static const G4float GREfactor = 5.0f;
  if(fExcEnergy >= (G4double)(GREfactor*GRWidth[A] + GREnergy[A])) { 
    return fProbability; 
  }
  // probability computed assuming continium transitions
  // VI: continium transition are limited only to final states
  //     below Fermi energy (this approach needs further evaluation)
  G4double emax = std::max(0.0, nucleus->ComputeGroundStateMass(Z, A-1) 
    + CLHEP::neutron_mass_c2 - nucleus->GetGroundStateMass());

  // max energy level for continues transition
  emax = std::min(emax, fExcEnergy);
  const G4double eexcfac = 0.99;
  if(0.0 == emax || fExcEnergy*eexcfac <= emax) { emax = fExcEnergy*eexcfac; }

  fStep = emax;
  const G4double MaxDeltaEnergy = CLHEP::MeV;
  fPoints = std::min((G4int)(fStep/MaxDeltaEnergy) + 2, MAXDEPOINT);
  fStep /= ((G4double)(fPoints - 1));
  if(fVerbose > 2) {
    G4cout << "Emax= " << emax << " Npoints= " << fPoints 
           << "  Eex= " << fExcEnergy << G4endl;
  }
  // integrate probabilities
  G4double eres = (G4double)GREnergy[A];
  G4double wres = (G4double)GRWidth[A];
  G4double eres2= eres*eres;
  G4double wres2= wres*wres;
  G4double levelDensity = fNuclearLevelData->GetLevelDensity(Z,A,fExcEnergy);
  G4double xsqr = std::sqrt(levelDensity*fExcEnergy);

  G4double egam    = fExcEnergy;
  G4double gammaE2 = egam*egam;
  G4double gammaR2 = gammaE2*wres2;
  G4double egdp2   = gammaE2 - eres2;

  G4double p0 = G4Exp(-2.0*xsqr)*gammaR2*gammaE2/(egdp2*egdp2 + gammaR2);
  G4double p1(0.0);

  for(G4int i=1; i<fPoints; ++i) {
    egam -= fStep;
    gammaE2 = egam*egam;
    gammaR2 = gammaE2*wres2;
    egdp2   = gammaE2 - eres2;
    p1 = G4Exp(2.0*(std::sqrt(levelDensity*std::abs(fExcEnergy - egam)) - xsqr))
      *gammaR2*gammaE2/(egdp2*egdp2 + gammaR2);
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
  fProbability *= fStep*NormC*A;
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
  if(fLevelManager) { 
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
  G4double eexc = nucleus->GetExcitationEnergy();
  if(eexc <= fTolerance) { return result; }

  InitialiseLevelManager(nucleus->GetZ_asInt(), nucleus->GetA_asInt());
  nucleus->SetLongLived(false);

  G4double time = nucleus->GetCreationTime();

  G4double efinal = 0.0;
  G4double ratio  = 0.0;
  vShellNumber    = -1;
  G4int  JP1      = 0;
  G4int  JP2      = 0;
  G4int  multiP   = 0;
  G4bool isGamma  = true;
  G4bool isDiscrete = false;

  const G4NucLevel* level = nullptr;
  std::size_t ntrans = 0;

  if(fVerbose > 2) {
    G4cout << "GenerateGamma: " << " Eex= " << eexc
           << " Eexmax= " << fLevelEnergyMax << G4endl;
  }
  // initial discrete state
  if(nullptr != fLevelManager && eexc <= fLevelEnergyMax + fTolerance) {
    fIndex = fLevelManager->NearestLevelIndex(eexc);
    G4double elevel = fLevelManager->LevelEnergy(fIndex); 
    isDiscrete = (std::abs(elevel - eexc) < fTolerance);
    if(fVerbose > 2) {
      G4cout << "              index= " << fIndex 
             << " lTime= " << fLevelManager->LifeTime(fIndex) << G4endl;
    }
    if(isDiscrete && 0 < fIndex) {
      // for discrete transition  
      level = fLevelManager->GetLevel(fIndex);
      if(nullptr != level) { 
        if(fVerbose > 2) {
          G4cout << "              ntrans= " << ntrans << " JP= " << JP1
                 << " RDM: " << fRDM << G4endl;
        }
        ntrans = level->NumberOfTransitions();
	// for floating level check levels with the same energy
        if(fLevelManager->FloatingLevel(fIndex) > 0 && 0 == ntrans &&
	   std::abs(elevel - fLevelManager->LevelEnergy(fIndex-1)) < fTolerance) {
	  auto newlevel = fLevelManager->GetLevel(fIndex-1);
	  if(nullptr != newlevel && newlevel->NumberOfTransitions() > 0) {
            --fIndex;
	    level = newlevel;
	    ntrans = level->NumberOfTransitions();
	  }
	}
        JP1 = fLevelManager->SpinTwo(fIndex); 
      }
    }
    // if a level has no defined transitions
    if(0 == ntrans) { isDiscrete = false; }
  }
  if(fVerbose > 2) {
    G4long prec = G4cout.precision(4);
    G4cout << "GenerateGamma: Z= " << nucleus->GetZ_asInt()
           << " A= " << nucleus->GetA_asInt() 
           << " Exc= " << eexc << " Emax= " 
           << fLevelEnergyMax << " idx= " << fIndex
           << " fCode= " << fCode << " fPoints= " << fPoints
           << " Ntr= " << ntrans << " discrete: " << isDiscrete
           << " fProb= " << fProbability << G4endl;
    G4cout.precision(prec);
  }

  // continues part
  if(!isDiscrete) {
    // we compare current excitation versus value used for probability 
    // computation and also Z and A used for probability computation 
    if(fCode != 1000*theZ + theA || eexc != fExcEnergy) { 
      GetEmissionProbability(nucleus); 
    }
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
    // final discrete level
    if(fVerbose > 2) {
      G4cout << "Continues proposes Efinal= " << efinal << G4endl;
    }
    if(nullptr != fLevelManager) {
      if(efinal < fLevelEnergyMax) {
        fIndex = fLevelManager->NearestLevelIndex(efinal, fIndex);
        efinal = fLevelManager->LevelEnergy(fIndex);
        // protection - take level below
        if(efinal >= eexc && 0 < fIndex) {
          --fIndex;
          efinal = fLevelManager->LevelEnergy(fIndex);
        } 
        nucleus->SetFloatingLevelNumber(fLevelManager->FloatingLevel(fIndex));

        // not allowed to have final energy above max energy
        // if G4LevelManager exist
      } else {
        efinal = fLevelEnergyMax;
        fIndex = fLevelManager->NearestLevelIndex(efinal, fIndex);
      }
    }
    if(fVerbose > 2) {
      G4cout << "Continues emission efinal(MeV)= " << efinal << G4endl; 
    }
    //discrete part ground state
  } else if(0 == fIndex) {
    G4bool isLL = false;
    if(nullptr != fLevelManager) {
      G4double ltime = fLevelManager->LifeTime(0);
      if(ltime < 0.0 || ltime > fMaxLifeTime) { isLL = true; }
    }
    nucleus->SetLongLived(isLL);
    return result;

    //discrete part
  } else {
 
    if(fVerbose > 2) {
      G4cout << "Discrete emission from level Index= " << fIndex 
             << " Elevel= " << fLevelManager->LevelEnergy(fIndex)
             << " Ltime= " << fLevelManager->LifeTime(fIndex)
             << " LtimeMax= " << fMaxLifeTime
             << "  RDM= " << fRDM << "  ICM= " << fICM << G4endl;
    }

    // stable fragment has life time -1 or above the limit
    // if is called from the radioactive decay the life time is not checked
    G4double ltime = fLevelManager->LifeTime(fIndex);
    if(ltime < 0.0 || (!fRDM && ltime > fMaxLifeTime)) {
      nucleus->SetLongLived(true);
      return result;
    }

    std::size_t idx = 0;
    if(1 < ntrans) {
      idx = level->SampleGammaTransition(G4UniformRand());
    }
    if(fVerbose > 2) {
      G4cout << "Ntrans= " << ntrans << " idx= " << idx
	     << " ICM= " << fICM << "  JP1= " << JP1 << G4endl;
    }
    G4double prob = level->GammaProbability(idx);
    // prob = 0 means that there is only internal conversion
    if(fICM && prob < 1.0) {
      G4double rndm = G4UniformRand();
      if(rndm > prob) {
	isGamma = false;
	rndm = (rndm - prob)/(1.0 - prob);
	vShellNumber = level->SampleShell(idx, rndm);
      }
    }
    // it is discrete transition with possible gamma correlation
    ratio  = level->MultipolarityRatio(idx);
    multiP = level->TransitionType(idx);
    fIndex = level->FinalExcitationIndex(idx);
    JP2 = fLevelManager->SpinTwo(fIndex); 

    // final energy and time
    efinal = fLevelManager->LevelEnergy(fIndex);
    // time is sampled if decay not prompt and this class called not 
    // from radioactive decay and isomer production is enabled 
    if(fSampleTime && ltime > 0.0) { 
      time -= ltime*G4Log(G4UniformRand()); 
    }
    nucleus->SetFloatingLevelNumber(fLevelManager->FloatingLevel(fIndex));
  }

  G4bool isLL = false;
  if(nullptr != fLevelManager) {
    G4double ltime = fLevelManager->LifeTime(fIndex);
    if(ltime < 0.0 || ltime > fMaxLifeTime) { isLL = true; }
  }
  nucleus->SetLongLived(isLL);

  // protection for floating levels
  if(std::abs(efinal - eexc) <= fTolerance) { return result; }

  result = fTransition->SampleTransition(nucleus, efinal, ratio, JP1,
                                         JP2, multiP, vShellNumber, 
                                         isDiscrete, isGamma);
  if(nullptr != result) { result->SetCreationTime(time); }

  // updated residual nucleus
  nucleus->SetCreationTime(time);
  nucleus->SetSpin(0.5*JP2);
  if(nullptr != fPolarization) { fPolarization->SetExcitationEnergy(efinal); }

  // ignore the floating levels with zero energy and create ground state
  if(efinal == 0.0 && fIndex > 0) {
    fIndex = 0;
    nucleus->SetFloatingLevelNumber(fLevelManager->FloatingLevel(0));
  }
      
  if(fVerbose > 2) { 
    G4cout << "Final level E= " << efinal << " time= " << time 
           << " idxFinal= " << fIndex << " isDiscrete: " << isDiscrete
           << " isGamma: " << isGamma << " multiP= " << multiP 
           << " shell= " << vShellNumber 
           << " JP1= " << JP1 << " JP2= " << JP2 << G4endl;
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
  
