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
//
// Geant4 header G4HadParticles
//
// Author V.Ivanchenko 09.05.2020
//
// Collections of hadrons pdg codes
//

#ifndef G4HadParticles_h
#define G4HadParticles_h 1

#include "globals.hh"
#include <vector>

class G4HadParticles
{
public:

  // several vectors of PDG codes for hadron physics
  static const std::vector<G4int>& GetLightHadrons();
  static const std::vector<G4int>& GetHyperons();
  static const std::vector<G4int>& GetAntiHyperons();
  static const std::vector<G4int>& GetKaons();
  static const std::vector<G4int>& GetBCHadrons();
  static const std::vector<G4int>& GetLightIons();
  static const std::vector<G4int>& GetLightAntiIons();

  // several vectors of PDG codes for EM physics
  static const std::vector<G4int>& GetHeavyChargedParticles();
  static const std::vector<G4int>& GetBCChargedHadrons();

private:

  static const std::vector<G4int> sLightHadrons;
  static const std::vector<G4int> sHyperons;
  static const std::vector<G4int> sAntiHyperons;
  static const std::vector<G4int> sKaons;
  static const std::vector<G4int> sBCHadrons;
  static const std::vector<G4int> sLightIons;
  static const std::vector<G4int> sLightAntiIons;
  static const std::vector<G4int> sHeavyChargedPart;
  static const std::vector<G4int> sBCChargedHadrons;
};

#endif


