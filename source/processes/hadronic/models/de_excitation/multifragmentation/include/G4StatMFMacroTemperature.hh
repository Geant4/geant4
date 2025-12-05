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

#ifndef G4StatMFMacroTemperature_h
#define G4StatMFMacroTemperature_h 1

#include <vector>
#include "globals.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4FunctionSolver.hh"

class G4StatMFMacroChemicalPotential;

class G4StatMFMacroTemperature {

public:

  G4StatMFMacroTemperature();
  ~G4StatMFMacroTemperature();

  void Initialise(const G4int anA, const G4int aZ, 
		  const G4double ExEnergy, const G4double FreeE0, 
		  const G4double kappa, 
		  std::vector<G4VStatMFMacroCluster*>* ClusterVector);
	
   
  G4double Function(G4double T)
  { return (fExEnergy - FragsExcitEnergy(T)); }	

  // copy constructor
  G4StatMFMacroTemperature(const G4StatMFMacroTemperature&) = delete;
  G4StatMFMacroTemperature& operator=
  (const G4StatMFMacroTemperature& right) = delete;
  G4bool operator==(const G4StatMFMacroTemperature& right) const = delete;
  G4bool operator!=(const G4StatMFMacroTemperature& right) const = delete;

  G4double GetMeanMultiplicity(void) const {return fMeanMultiplicity;}
	
  G4double GetChemicalPotentialMu(void) const {return fChemPotentialMu;}

  G4double GetChemicalPotentialNu(void) const {return fChemPotentialNu;}

  G4double GetTemperature(void) const {return fMeanTemperature;}

  G4double GetEntropy(void) const {return fMeanEntropy;}

  G4double CalcTemperature(void);

private:
	
  G4double FragsExcitEnergy(const G4double T);

  void CalcChemicalPotentialNu(const G4double T);

  G4FunctionSolver<G4StatMFMacroTemperature>* fSolver;
  G4StatMFMacroChemicalPotential* theChemPot;

  G4int theA{0};
  G4int theZ{0};
  G4double fExEnergy{0.0};
  G4double fFreeInternalE0{0.0};
  G4double fKappa{0.0};
  G4double fMeanMultiplicity{0.0};
  G4double fMeanTemperature{0.0};
  G4double fChemPotentialMu{0.0};
  G4double fChemPotentialNu{0.0};
  G4double fMeanEntropy{0.0};
	
  std::vector<G4VStatMFMacroCluster*>* fClusters{nullptr};
};

#endif
