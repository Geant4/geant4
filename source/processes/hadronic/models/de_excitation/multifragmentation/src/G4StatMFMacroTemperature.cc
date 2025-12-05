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
// Modified:
// 25.07.08 I.Pshenichnov (in collaboration with Alexander Botvina and Igor 
//          Mishustin (FIAS, Frankfurt, INR, Moscow and Kurchatov Institute, 
//          Moscow, pshenich@fias.uni-frankfurt.de) make algorithm closer to
//          original MF model
// 16.04.10 V.Ivanchenko improved logic of solving equation for temperature
//          to protect code from rare unwanted exception; moved constructor 
//          and destructor to source  
// 28.10.10 V.Ivanchenko defined members in constructor and cleaned up
// 13.08.2025 V.Ivanchenko rewrite

#include "G4StatMFMacroTemperature.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFMacroChemicalPotential.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"

namespace {
  const G4double t1 = 1*CLHEP::MeV;
  const G4double t2 = 50*CLHEP::MeV;
}

G4StatMFMacroTemperature::G4StatMFMacroTemperature()
{
  fSolver = new G4FunctionSolver<G4StatMFMacroTemperature>(this, 100, 5.e-4);
  fSolver->SetIntervalLimits(t1, t2);
  theChemPot = new G4StatMFMacroChemicalPotential();
}

G4StatMFMacroTemperature::~G4StatMFMacroTemperature()
{
  delete fSolver;
  delete theChemPot;
}

void G4StatMFMacroTemperature::Initialise(const G4int anA, const G4int aZ, 
					  const G4double ExEnergy,
					  const G4double FreeE0,
					  const G4double kappa, 
					  std::vector<G4VStatMFMacroCluster*>* v)
{
  theA = anA;
  theZ = aZ;
  fExEnergy = ExEnergy;
  fFreeInternalE0 = FreeE0;
  fKappa = kappa;
  fClusters = v;
}

G4double G4StatMFMacroTemperature::CalcTemperature(void) 
{
  fMeanTemperature = std::max(std::min(std::sqrt(fExEnergy/(theA*0.12)), t2), t1);
  fSolver->FindRoot(fMeanTemperature);
  return fMeanTemperature;
}

G4double G4StatMFMacroTemperature::FragsExcitEnergy(const G4double T)
// Calculates excitation energy per nucleon and summed fragment 
// multiplicity and entropy
{
  // Model Parameters
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double R0 = G4StatMFParameters::Getr0()*g4calc->Z13(theA);
  G4double R = R0*g4calc->A13(1.0 + G4StatMFParameters::GetKappaCoulomb());
  G4double FreeVol = fKappa*(4.*CLHEP::pi/3.)*R0*R0*R0; 
 
  // Calculate Chemical potentials
  CalcChemicalPotentialNu(T);

  // Average total fragment energy and mean entropy
  G4double AverageEnergy = 0.0;
  fMeanEntropy = 0.0;
  for (auto const ptr : *fClusters) {
    fMeanEntropy += ptr->CalcEntropy(T, FreeVol);
    AverageEnergy += ptr->GetMeanMultiplicity() * ptr->CalcEnergy(T);
  }
    
  // Add Coulomb energy			
  AverageEnergy += 0.6*CLHEP::elm_coupling*(theZ*theZ)/R;

  // Excitation energy per nucleon
  return AverageEnergy - fFreeInternalE0;
}

void G4StatMFMacroTemperature::CalcChemicalPotentialNu(const G4double T)
// Calculates the chemical potential \nu 
{
  theChemPot->Initialise(theA, theZ, fKappa, T, fClusters);

  fChemPotentialNu = theChemPot->CalcChemicalPotentialNu();
  fChemPotentialMu = theChemPot->GetChemicalPotentialMu();
  fMeanMultiplicity = theChemPot->GetMeanMultiplicity();		
}


