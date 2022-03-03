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
/// \file Event.hh
/// \brief Definition of the CaTS::Event class

#pragma once

#include <G4String.hh>
#include <vector>
#include <map>
#include "G4Types.hh"
class G4VHit;

class Event
{
 public:
  Event()          = default;
  virtual ~Event() = default;
  inline void SetEventNr(G4int i) { fEvtNum = i; }
  inline G4int GetEventNumber() const { return fEvtNum; }
  inline std::map<G4String, std::vector<G4VHit*>>* GetHCMap() { return &hcmap; }
  void Reset();

 private:
  G4int fEvtNum{ 0 };
  std::map<G4String, std::vector<G4VHit*>> hcmap;  // map of Hit Collections
};
