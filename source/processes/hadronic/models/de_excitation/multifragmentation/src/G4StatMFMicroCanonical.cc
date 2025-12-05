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
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modification: 13.08.2025 V.Ivanchenko rewrite

#include "G4StatMFMicroCanonical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"
#include "G4Pow.hh"

namespace
{
  constexpr G4int fMaxMultiplicity = 4;
  constexpr G4double t1 = 1*CLHEP::MeV;
  constexpr G4double t2 = 50*CLHEP::MeV;
}
  
// constructor
G4StatMFMicroCanonical::G4StatMFMicroCanonical() 
{
  fSolver = new G4FunctionSolver<G4StatMFMicroCanonical>(this, 100, 5.e-4);
  fSolver->SetIntervalLimits(t1, t2);
  fPartitionManagerVector.reserve(fMaxMultiplicity);
  g4calc = G4Pow::GetInstance();
}

// destructor
G4StatMFMicroCanonical::~G4StatMFMicroCanonical() 
{
  delete fSolver;
  if (!fPartitionManagerVector.empty()) {
    for (auto const & p : fPartitionManagerVector) { delete p; }
  }
}

void G4StatMFMicroCanonical::Initialise(const G4Fragment& theFragment) 
{
  fPartitionManagerVector.clear();
  // Excitation Energy 
  fExEnergy = theFragment.GetExcitationEnergy();

  A = theFragment.GetA_asInt();
  Z = theFragment.GetZ_asInt();
  A13 = g4calc->Z13(A);
  
  fInvLevelDensity = G4StatMFParameters::GetEpsilon0()*(1.0 + 3.0/G4double(A-1));
  
  fSymmetryTerm = G4StatMFParameters::GetGamma0()*(A - 2*Z)*(A - 2*Z)/(G4double)A;

  fCoulombTerm = elm_coupling*0.6*Z*Z/(G4StatMFParameters::Getr0()*A13);

  // Configuration temperature
  G4double TConf = std::sqrt(8.0*fExEnergy/(G4double)A);
  
  // Free internal energy at Temperature T = 0 (SurfaceTerm at T = 0)
  pFreeInternalE0 = -G4StatMFParameters::GetE0()*A + fSymmetryTerm  
    + G4StatMFParameters::GetBeta0()*A13*A13 + fCoulombTerm;
  
  //G4cout << "Tconf=" <<  TConf << " freeE=" << pFreeInternalE0 << G4endl;
    
  // Mean breakup multiplicity
  pMeanMultiplicity = 0.0;
  
  // Mean channel temperature
  pMeanTemperature = 0.0;
  
  // Mean channel entropy
  pMeanEntropy = 0.0;
  
  // Calculate entropy of compound nucleus
  G4double SCompoundNucleus = CalcEntropyOfCompoundNucleus(TConf);
  
  // Statistical weight of compound nucleus
  fWCompoundNucleus = 1.0; 
  
  // Statistical weight
  G4double W = fWCompoundNucleus;
  // Maximal fragment multiplicity allowed in direct simulation  

  for (G4int im = 2; im <= fMaxMultiplicity; ++im) {
    auto ptr = new G4StatMFMicroManager(theFragment, im, pFreeInternalE0, SCompoundNucleus);
    fPartitionManagerVector.push_back(ptr);
    W += ptr->GetProbability();
  }
  
  // Normalization of statistical weights
  for (auto & ptr : fPartitionManagerVector) {
    ptr->Normalize(W);
    pMeanMultiplicity += ptr->GetMeanMultiplicity();
    pMeanTemperature += ptr->GetMeanTemperature();
    pMeanEntropy += ptr->GetMeanEntropy();
  }

  fWCompoundNucleus /= W;
  
  pMeanMultiplicity += fWCompoundNucleus;
  pMeanTemperature += TConf * fWCompoundNucleus;
  pMeanEntropy += SCompoundNucleus * fWCompoundNucleus;
}

G4double G4StatMFMicroCanonical::CalcFreeInternalEnergy(G4double T)
{
  G4double VolumeTerm = (-G4StatMFParameters::GetE0()+T*T/fInvLevelDensity)*A;
  G4double SurfaceTerm = (G4StatMFParameters::Beta(T) - T*G4StatMFParameters::DBetaDT(T))*A13*A13;
  G4double sum = VolumeTerm + fSymmetryTerm + SurfaceTerm + fCoulombTerm;
  
  // G4cout << "G4StatMFMicroCanonical::CalcFreeInternalEnergy " << sum
  // 	 << " " << VolumeTerm << " " << fSymmetryTerm << " " << SurfaceTerm
  //	 <<  " " << fCoulombTerm << G4endl;
  return sum;
}

G4double G4StatMFMicroCanonical::CalcEntropyOfCompoundNucleus(G4double& TConf)
  // Calculates Temperature and Entropy of compound nucleus
{
  G4double T = std::max(std::min(std::max(TConf,std::sqrt(fExEnergy/(A*0.125))), t2), t1);\
  fSolver->FindRoot(T);
  TConf = T;
  // G4cout << "=== FindRoot T= " << T << G4endl;
  auto S = (2*A)*T/fInvLevelDensity - G4StatMFParameters::DBetaDT(T)*A13*A13;
  return S;
}

G4StatMFChannel*  G4StatMFMicroCanonical::ChooseAandZ(const G4Fragment& theFragment)
{
  // Choice of fragment atomic numbers and charges 
  // We choose a multiplicity (1,2,3,...) and then a channel
  G4int AA = theFragment.GetA_asInt();
  G4int ZZ = theFragment.GetZ_asInt();
  
  if (G4UniformRand() < fWCompoundNucleus) { 
	
    G4StatMFChannel * aChannel = new G4StatMFChannel;
    aChannel->CreateFragment(AA, ZZ);
    return aChannel;
	
  } else {
    G4double rand = G4UniformRand();
    G4double AccumWeight = fWCompoundNucleus;
    for (auto & ptr : fPartitionManagerVector) {
      AccumWeight += ptr->GetProbability();
      if (rand <= AccumWeight) {
	return ptr->ChooseChannel(A, Z, pMeanTemperature);
      }
    }
  }
  return nullptr;
}
