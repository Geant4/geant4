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

#ifndef G4StatMFMacroChemicalPotential_h
#define G4StatMFMacroChemicalPotential_h 1

#include <vector>
#include "G4VStatMFMacroCluster.hh"
#include "G4FunctionSolver.hh"

class G4StatMFMacroMultiplicity;

class G4StatMFMacroChemicalPotential {

public:

  G4StatMFMacroChemicalPotential();

  ~G4StatMFMacroChemicalPotential();

  void Initialise(const G4int anA, const G4int aZ,
		  const G4double kappa, const G4double temp, 
		  std::vector<G4VStatMFMacroCluster*>* cVector);

  G4double Function(G4double nu)
  { return (theZ - CalcMeanZ(nu)); }

  G4double CalcChemicalPotentialNu();

  G4double GetMeanMultiplicity() const {return fMeanMultiplicity;}	
  G4double GetChemicalPotentialMu() const {return fChemPotentialMu;}
  G4double GetChemicalPotentialNu() const {return fChemPotentialNu;}

  G4StatMFMacroChemicalPotential(const G4StatMFMacroChemicalPotential &) = delete;
  G4StatMFMacroChemicalPotential& operator=
  (const G4StatMFMacroChemicalPotential & right) = delete;
  G4bool operator==(const G4StatMFMacroChemicalPotential & right) const = delete;
  G4bool operator!=(const G4StatMFMacroChemicalPotential & right) const = delete;

private:
	
  G4double CalcMeanZ(const G4double nu);
  void CalcChemicalPotentialMu(const G4double nu);

  G4FunctionSolver<G4StatMFMacroChemicalPotential>* fSolver;
  G4StatMFMacroMultiplicity* theMultip;
  
  G4int theA{0};
  G4int theZ{0};
  G4double fKappa{0.0};
  G4double fMeanTemperature{0.0};
  G4double fMeanMultiplicity{0.0};
  G4double fChemPotentialMu{0.0};
  G4double fChemPotentialNu{0.0};
	
  std::vector<G4VStatMFMacroCluster*>* fClusters{nullptr}; 
};
#endif
