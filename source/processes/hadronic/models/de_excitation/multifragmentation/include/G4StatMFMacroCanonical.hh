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

#ifndef G4StatMFMacroCanonical_h
#define G4StatMFMacroCanonical_h 1

#include "G4Fragment.hh"
#include "G4StatMFFragment.hh"
#include "G4VStatMFEnsemble.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4StatMFMacroNucleon.hh"
#include "G4StatMFMacroBiNucleon.hh"
#include "G4StatMFMacroTriNucleon.hh"
#include "G4StatMFMacroTetraNucleon.hh"
#include "G4StatMFMacroMultiNucleon.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"
#include "G4StatMFMacroTemperature.hh"
#include "Randomize.hh"


class G4StatMFMacroCanonical : public G4VStatMFEnsemble {

public:

  G4StatMFMacroCanonical();

  ~G4StatMFMacroCanonical() override;

  // Initialise for a given G4Fragment
  void Initialise(const G4Fragment& theFragment) override;

  // Choice of the channel
  G4StatMFChannel* ChooseAandZ(const G4Fragment &theFragment) override;

  G4StatMFMacroCanonical(const G4StatMFMacroCanonical&) = delete;
  G4StatMFMacroCanonical& operator=(const G4StatMFMacroCanonical& right) = delete;
  G4bool operator==(const G4StatMFMacroCanonical& right) const = delete;
  G4bool operator!=(const G4StatMFMacroCanonical& right) const = delete;

private:

  // Determines fragments multiplicities and compute total fragment multiplicity
  G4double ChooseA(G4int A, std::vector<G4int>& ANumbers);
	
  // Samples charges of fragments
  G4StatMFChannel* ChooseZ(G4int Z, std::vector<G4int>& FragmentsA);

  G4StatMFMacroTemperature* theTemp{nullptr};
  
  // Chemical Potential \mu
  G4double fChemPotentialMu{0.0};

  // Chemical Potential \nu
  G4double fChemPotentialNu{0.0};

  // Parameter Kappa
  G4double fKappa{0.0};

  // Clusters
  std::vector<G4VStatMFMacroCluster*> fClusters;
  std::vector<G4double> fAcumMultiplicity;
};

#endif
