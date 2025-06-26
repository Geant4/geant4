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
// GEANT4 Class file
//
// File name: G4PreCompoundInterface
//
// Author:  V.Ivantchenko, 20 January 2025
//

#include "G4PreCompoundInterface.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PreCompoundEmissionInt.hh"
#include "G4PreCompoundTransitionsInt.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4NucleiProperties.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4LorentzVector.hh"
#include "G4Exp.hh"
#include "G4PhysicsModelCatalog.hh"

////////////////////////////////////////////////////////////////////////////////

G4PreCompoundInterface::G4PreCompoundInterface() 
  : G4VPreCompoundModel(new G4ExcitationHandler(),"PRECO_I")
{
  fNuclData = G4NuclearLevelData::GetInstance();
}

////////////////////////////////////////////////////////////////////////////////

G4PreCompoundInterface::~G4PreCompoundInterface() 
{
  delete theEmission;
  delete theTransition;
  delete GetExcitationHandler();
}

////////////////////////////////////////////////////////////////////////////////

void G4PreCompoundInterface::BuildPhysicsTable(const G4ParticleDefinition&) 
{
  InitialiseModel();
}

////////////////////////////////////////////////////////////////////////////////

void G4PreCompoundInterface::InitialiseModel() 
{
  if (isInitialised) { return; }
  isInitialised = true;

  G4DeexPrecoParameters* param = fNuclData->GetParameters();

  fLowLimitExc = param->GetPrecoLowEnergy();
  fHighLimitExc = param->GetPrecoHighEnergy();
  fVerbose = param->GetVerbose();

  if (param->PrecoDummy()) {
    isActive = false;
  } else {
    theEmission = new G4PreCompoundEmissionInt(fVerbose);
    theEmission->SetOPTxs(param->GetPrecoModelType());
    theTransition = new G4PreCompoundTransitionsInt(fVerbose);
  }

  GetExcitationHandler()->Initialise();
}

////////////////////////////////////////////////////////////////////////////////

G4ReactionProductVector* G4PreCompoundInterface::DeExcite(G4Fragment& frag)
{
  G4ReactionProductVector* res = new G4ReactionProductVector;
  if (!isInitialised) { InitialiseModel(); }

  // decays by de-excitation
  G4double U = frag.GetExcitationEnergy();
  G4int Z = frag.GetZ_asInt(); 
  G4int A = frag.GetA_asInt();
  if (1 < fVerbose) {
    G4cout << "### G4PreCompoundInterface::DeExcite Z=" << Z << " A=" << A
           << " U(MeV)=" << U << G4endl;
  }
  if (!isActive || Z < minZ || A < minA || 
      U < fLowLimitExc*A || U > A*fHighLimitExc ||
      0 < frag.GetNumberOfLambdas()) {
    PerformEquilibriumEmission(frag, res);
    return res;
  }
  // decays by precompound model 
  BreakUpFragment(frag, res);
  return res;
}

void G4PreCompoundInterface::BreakUpFragment(G4Fragment& frag,
                                             G4ReactionProductVector* res)
{
  // check initial fragment and number of excitons is not defined
  G4double U = frag.GetExcitationEnergy();
  G4int Z = frag.GetZ_asInt(); 
  G4int A = frag.GetA_asInt();
  if (Z < minZ || A < minA || U < fLowLimitExc*A) {
    PerformEquilibriumEmission(frag, res);
    return;
  }

  const G4double ldfact = 3.0/CLHEP::pi2;
  const G4double eperex = 20.0*CLHEP::MeV;
  
  // number of excitons may be defined or not defined
  G4int np = frag.GetNumberOfParticles();
  G4int nz = frag.GetNumberOfCharged();
  if (0 == np) {
    np = G4lrint(U/eperex);
    np = std::min(np, A);
    frag.SetNumberOfParticles(np);
    nz = G4lrint((np*Z)/(G4double)A);
    nz = std::min(std::min(nz, np), Z);
    frag.SetNumberOfExcitedParticle(np, nz);
  }

  // main loop over fragments
  for (G4int i=0; i<50; ++i) {
    if (Z < minZ || A < minA || U < fLowLimitExc*A) {
      break;
    }

    // eqNum is the number of particle, not excitons
    G4int eqNum =
      G4lrint(std::sqrt(ldfact*U*fNuclData->GetLevelDensity(Z, A, U)));        
    if (2 < fVerbose) {
      G4cout << "   1st loop " << i << ". Z=" << Z << " A=" << A
             << " U(MeV)=" << U << " Npart=" << np << " Nch=" << nz
	     << " eqExcitationNumber=" << eqNum << G4endl;
    }
    if (np <= eqNum) { break; }

    // Loop for transitions, it is performed while there are 
    // preequilibrium transitions and is completed by emission
    G4bool isTransition = false;
    G4bool isEquilibrium = false;
    
    for (G4int j=0; j<20; ++j) {
      G4double transProbability =
	theTransition->CalculateProbability(frag);
      G4double P1 = theTransition->GetTransitionProb1();
      G4double P2 = theTransition->GetTransitionProb2();
      G4double P3 = theTransition->GetTransitionProb3();
      if (2 < fVerbose) {
        G4cout << "   2nd loop " << j << ". Npart=" << np << " P1=" << P1
             << " P2=" << P2 << " P3=" << P3 << G4endl;
      }
      if (np <= eqNum || P1 <= P2+P3) {
	isEquilibrium = true;
	break;
      }
      
      G4double emissionProbability =
	theEmission->GetTotalProbability(frag);

      // Sum of all probabilities
      G4double totalProb = emissionProbability + transProbability;
            
      // Select subprocess
      if (totalProb*G4UniformRand() > emissionProbability) {
	isTransition = true;
	theTransition->PerformTransition(frag);
	isEquilibrium = (np <= eqNum);
      } else {
	isTransition = false;
	auto product = theEmission->PerformEmission(frag);
	res->push_back(product);

	// new parameters of the residual fragment
	U = frag.GetExcitationEnergy();
	Z = frag.GetZ_asInt(); 
	A = frag.GetA_asInt();
      }
      np = frag.GetNumberOfParticles();
      nz = frag.GetNumberOfCharged();
      if (!isTransition || isEquilibrium) { break; }
    }
    if (isEquilibrium) { break; }
  }
  PerformEquilibriumEmission(frag, res);
}

////////////////////////////////////////////////////////////////////////////////
//       Documentation
////////////////////////////////////////////////////////////////////////////////

void G4PreCompoundInterface::ModelDescription(std::ostream& outFile) const
{
  outFile 
    << "The GEANT4 precompound model is considered as an extension of the\n"
    <<	"hadron kinetic model. It gives a possibility to extend the low energy range\n"
    <<	"of the hadron kinetic model for nucleon-nucleus inelastic collision and it \n"
    <<	"provides a ”smooth” transition from kinetic stage of reaction described by the\n"
    <<	"hadron kinetic model to the equilibrium stage of reaction described by the\n"
    <<	"equilibrium deexcitation models.\n"
    <<	"The initial information for calculation of pre-compound nuclear stage\n"
    <<	"consists of the atomic mass number A, charge Z of residual nucleus, its\n"
    <<	"four momentum P0 , excitation energy U and number of excitons n, which equals\n"
    <<	"the sum of the number of particles p (from them p_Z are charged) and the number of\n"
    <<	"holes h.\n"
    <<	"At the preequilibrium stage of reaction, we follow the exciton model approach in ref. [1],\n"
    <<	"taking into account the competition among all possible nuclear transitions\n"
    <<	"with ∆n = +2, −2, 0 (which are defined by their associated transition probabilities) and\n"
    <<	"the emission of neutrons, protons, deuterons, thritium and helium nuclei (also defined by\n"
    <<	"their associated emission  probabilities according to exciton model)\n"
    <<	"\n"
    <<	"[1] K.K. Gudima, S.G. Mashnik, V.D. Toneev, Nucl. Phys. A401 329 (1983)\n"
    <<  "\n";
}

void G4PreCompoundInterface::DeExciteModelDescription(std::ostream& outFile) const
{
  outFile << "description of precompound model as used with DeExcite()" << "\n";
}
  
////////////////////////////////////////////////////////////////////////////////
