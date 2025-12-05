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
//          Moscow, pshenich@fias.uni-frankfurt.de) additional checks in
//          solver of equation for the chemical potential
//
// 13.08.2025 V.Ivanchenko rewrite

#include "G4StatMFMacroMultiplicity.hh"
#include "G4PhysicalConstants.hh"
#include "G4StatMFParameters.hh"
#include "G4Pow.hh"

G4StatMFMacroMultiplicity::G4StatMFMacroMultiplicity()
{
  fSolver = new G4FunctionSolver<G4StatMFMacroMultiplicity>(this, 100, 5.e-4);
}

G4StatMFMacroMultiplicity::~G4StatMFMacroMultiplicity()
{
  delete fSolver;
}

void G4StatMFMacroMultiplicity::Initialise(const G4int anA, 
					   const G4double kappa, 
					   const G4double temp, 
					   const G4double nu,
					   std::vector<G4VStatMFMacroCluster*>* v)
{
  A = anA;
  theA = anA;
  fKappa = kappa;
  fMeanTemperature = temp;
  fChemPotentialNu = nu;
  fClusters = v;
}

G4double G4StatMFMacroMultiplicity::CalcChemicalPotentialMu() 
    // Calculate Chemical potential \mu
    // For that is necesary to calculate mean multiplicities
{
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double CP = G4StatMFParameters::GetCoulomb();

  // starting value for chemical potential \mu
  // it is the derivative of F(T,V)-\nu*Z w.r.t. Af in Af=5
  G4double ZA5 = (*fClusters)[4]->GetZARatio();
  G4double ILD5 = (*fClusters)[4]->GetInvLevelDensity();
  fChemPotentialMu = -G4StatMFParameters::GetE0()-
    fMeanTemperature*fMeanTemperature/ILD5 -
    fChemPotentialNu*ZA5 + 
    G4StatMFParameters::GetGamma0()*(1.0-2.0*ZA5)*(1.0-2.0*ZA5) +
    (2.0/3.0)*G4StatMFParameters::Beta(fMeanTemperature)/g4calc->Z13(5) +
    (5.0/3.0)*CP*ZA5*ZA5*g4calc->Z23(5) -
    1.5*fMeanTemperature/5.0;
		
  G4double ChemPa = fChemPotentialMu;
  if (ChemPa > 10*fMeanTemperature) { ChemPa = 10*fMeanTemperature; }
  G4double ChemPb = ChemPa - 0.5*std::abs(ChemPa);
  fSolver->SetIntervalLimits(ChemPa, ChemPb);
  fSolver->FindRoot(fChemPotentialMu); 
  return fChemPotentialMu;
}

G4double G4StatMFMacroMultiplicity::CalcMeanA(const G4double mu)
{
  G4double r0 = G4StatMFParameters::Getr0(); 
  G4double V0 = (4.0/3.0)*pi*theA*r0*r0*r0;

  G4double MeanA = 0.0;
	
  fMeanMultiplicity = 0.0;
	
  G4int n = 1;
  G4int nn = (G4int)fClusters->size();
  nn = std::min(nn, A);
  for (G4int i=0; i<nn; ++i) {
    G4double multip =
      (*fClusters)[i]->CalcMeanMultiplicity(V0*fKappa,mu,fChemPotentialNu,
					    fMeanTemperature);
    MeanA += multip*(++n);
    fMeanMultiplicity += multip;
  }

  return MeanA;
}
