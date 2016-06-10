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
// $Id: MCTruthConfig.hh 73446 2013-08-27 11:32:59Z gcosmo $
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

  void SetMinE(double e) {minE = e;}
  G4double GetMinE() const {return minE;}

  void SetParticleTypes(std::vector<G4int>& types) {particleTypes = types;}
  void AddParticleType(G4int type) {particleTypes.push_back(type);}
  std::vector<G4int>& GetParticleTypes() {return particleTypes;}

  void SetCreatorProcesses(std::vector<G4String>& processes)
       {creatorProcesses = processes;}
  void AddCreatorProcess(G4String& process)
       {creatorProcesses.push_back(process);}
  std::vector<G4String>& GetCreatorProcesses()
       {return creatorProcesses;}

private:

  G4double minE;
  std::vector<G4int> particleTypes;
  std::vector<G4String> creatorProcesses;

};

#endif // INCLUDE_MCTRUTHCONFIG_H
