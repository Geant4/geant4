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

// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the CaTS::PrimaryGeneratorAction class

#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include <G4String.hh>
#include <G4ios.hh>
#include <map>
#include <utility>
class G4GenericMessenger;
class G4VPrimaryGenerator;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();
  void GeneratePrimaries(G4Event*) final;
  inline void SetGenerator(G4String genname)
  {
    std::map<G4String, G4VPrimaryGenerator*>::iterator pos =
      gentypeMap.find(genname);
    if(pos != gentypeMap.end())
    {
      fcurrentGenerator     = pos->second;
      fcurrentGeneratorName = genname;
    }
    else
    {
      G4cout << "Primary  Generator: " << genname << " is not a valid Option"
             << G4endl;
      G4cout << "Valid Options are: ";
#if(G4VERSION_NUMBER > 1072)
      for(auto [genType, generator] : gentypeMap)
      {
        G4cout << genType << " ";
      }
      G4cout << G4endl;
#else
      for(std::map<G4String, G4VPrimaryGenerator*>::iterator ii =
            gentypeMap.begin();
          ii != gentypeMap.end(); ++ii)
      {
        G4cout << (*ii).first << " ";
      }
      G4cout << G4endl;
#endif
    }
  }
  inline G4VPrimaryGenerator* GetGenerator() const { return fcurrentGenerator; }
  inline G4String GetGeneratorName() const { return fcurrentGeneratorName; }
  inline void Print()
  {
    G4cout << "===================================================" << G4endl;
    G4cout << " Primary  Generator: " << fcurrentGeneratorName << G4endl;
    G4cout << "===================================================" << G4endl;
  }

 private:
  /// Define UI commands: choice of primary generators
  void DefineCommands();
  G4VPrimaryGenerator* fcurrentGenerator = nullptr;
  G4String fcurrentGeneratorName{ "particleGun" };
  PrimaryGeneratorAction& operator=(const PrimaryGeneratorAction& right);
  PrimaryGeneratorAction(const PrimaryGeneratorAction&);
  std::map<G4String, G4VPrimaryGenerator*> gentypeMap;
  /// Pointer to the messenger for UI commands
  G4GenericMessenger* fMessenger = nullptr;
};
