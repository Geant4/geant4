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

#include "G4StatMFMacroChemicalPotential.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFMacroMultiplicity.hh"
#include "G4PhysicalConstants.hh"
#include "G4Pow.hh"

G4StatMFMacroChemicalPotential::G4StatMFMacroChemicalPotential()
{
  theMultip = new G4StatMFMacroMultiplicity();
  fSolver = new G4FunctionSolver<G4StatMFMacroChemicalPotential>(this, 100, 5.e-4);
}


G4StatMFMacroChemicalPotential::~G4StatMFMacroChemicalPotential()
{
  delete fSolver;
  delete theMultip;
}

void G4StatMFMacroChemicalPotential::Initialise(
				     const G4int anA, const G4int aZ,
				     const G4double kappa,
				     const G4double temp, 
				     std::vector<G4VStatMFMacroCluster*>* v)
{
  theA = anA;
  theZ = aZ;
  fKappa = kappa;
  fMeanTemperature = temp;
  fClusters = v;
}

G4double G4StatMFMacroChemicalPotential::CalcChemicalPotentialNu()
//	Calculate Chemical potential \nu
{
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double CP = G4StatMFParameters::GetCoulomb();

  // Initial value for fChemPotentialNu	
  fChemPotentialNu = (theZ/(G4double)theA)*
    (8.0*G4StatMFParameters::GetGamma0() + 2.0*CP*g4calc->Z23(theA))
    - 4.0*G4StatMFParameters::GetGamma0();
		
  fSolver->SetIntervalLimits(0.5*fChemPotentialNu, 2*fChemPotentialNu);
  fSolver->FindRoot(fChemPotentialNu);
  return fChemPotentialNu;
}

G4double G4StatMFMacroChemicalPotential::CalcMeanZ(const G4double nu)
{
  CalcChemicalPotentialMu(nu);
  // This is important, the Z over A ratio for proton and neutron depends on the 
  // chemical potential Mu, while for the first guess for Chemical potential mu 
  // some values of Z over A ratio. This is the reason for that.
  
  G4double MeanZ = 0.0;
  G4int n = 0;
  G4int nn = (G4int)fClusters->size();
  nn = std::min(nn, theA);
  for (G4int i = 0; i < nn; ++i) {
    G4double x = (*fClusters)[i]->CalcZARatio(nu);
    MeanZ += (n++) * x * (*fClusters)[i]->GetMeanMultiplicity(); 
  }
  return MeanZ;
}

void G4StatMFMacroChemicalPotential::CalcChemicalPotentialMu(const G4double nu)
//	Calculate Chemical potential \mu
// For that is necesary to calculate mean multiplicities
{
  theMultip->Initialise(theA, fKappa, fMeanTemperature, nu, fClusters);
  fChemPotentialMu = theMultip->CalcChemicalPotentialMu();
  fMeanMultiplicity = theMultip->GetMeanMultiplicity();
}
