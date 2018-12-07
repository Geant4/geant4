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
/// \file eventgenerator/HepMC/MCTruth/include/MCTruthConfig.hh
/// \brief Definition of the MCTruthConfig class
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - MCTruthConfig class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
// Date  : 2006-03-06
//
// --------------------------------------------------------------
#ifndef INCLUDE_MCTRUTHCONFIG_H 
#define INCLUDE_MCTRUTHCONFIG_H 1

#include<vector>
#include<iostream>

#include "globals.hh"

class MCTruthConfig 
{
public: 

  MCTruthConfig(); 

  virtual ~MCTruthConfig();

  void SetMinE(double e) {fMinE = e;}
  G4double GetMinE() const {return fMinE;}

  void SetParticleTypes(std::vector<G4int>& types) {fParticleTypes = types;}
  void AddParticleType(G4int type) {fParticleTypes.push_back(type);}
  std::vector<G4int>& GetParticleTypes() {return fParticleTypes;}

  void SetCreatorProcesses(std::vector<G4String>& processes)
       {fCreatorProcesses = processes;}
  void AddCreatorProcess(G4String& process)
       {fCreatorProcesses.push_back(process);}
  std::vector<G4String>& GetCreatorProcesses()
       {return fCreatorProcesses;}

private:

  G4double fMinE;
  std::vector<G4int> fParticleTypes;
  std::vector<G4String> fCreatorProcesses;

};

#endif // INCLUDE_MCTRUTHCONFIG_H
