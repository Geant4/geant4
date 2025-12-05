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

#ifndef G4StatMFMacroMultiplicity_h
#define G4StatMFMacroMultiplicity_h 1

#include <vector>
#include "globals.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4FunctionSolver.hh"

class G4StatMFMacroMultiplicity {

public:

  G4StatMFMacroMultiplicity();

  ~G4StatMFMacroMultiplicity();

  void Initialise(const G4int anA, const G4double kappa, 
		  const G4double temp, const G4double nu,
		  std::vector<G4VStatMFMacroCluster*>* cVector);

  G4double Function(G4double mu)
  { return (theA - CalcMeanA(mu)); };

  G4double CalcChemicalPotentialMu();

  G4double GetMeanMultiplicity() const { return fMeanMultiplicity; }

  G4double GetChemicalPotentialMu() const { return fChemPotentialMu; }

  G4StatMFMacroMultiplicity(const G4StatMFMacroMultiplicity&) = delete;
  G4StatMFMacroMultiplicity& operator=
  (const G4StatMFMacroMultiplicity& right) = delete;
  G4bool operator==(const G4StatMFMacroMultiplicity& right) const = delete;
  G4bool operator!=(const G4StatMFMacroMultiplicity& right) const = delete;

private:
	
  G4double CalcMeanA(const G4double mu);

  G4int A{0};
  G4double theA{0};
  G4double fKappa{0.0};
  G4double fMeanTemperature{0.0};
  G4double fChemPotentialNu{0.0};

  G4double fMeanMultiplicity{0.0};
  G4double fChemPotentialMu{0.0};

  std::vector<G4VStatMFMacroCluster*>* fClusters{nullptr}; 
  G4FunctionSolver<G4StatMFMacroMultiplicity>* fSolver;
};

#endif
