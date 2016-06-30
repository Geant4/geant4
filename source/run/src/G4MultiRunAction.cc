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
// $Id: G4MultiRunAction.hh 90212 2016-01-27 18:33:12Z adotti $
//
//---------------------------------------------------------------
//
// G4MultiRunAction.hh
//
//   Created on: Jan 17, 2016
//       Author: adotti
//
// ---------------------------------------------------------------
//

#include "G4MultiRunAction.hh"
#include "G4Run.hh"
#include <algorithm>

G4Run* G4MultiRunAction::GenerateRun() {
  G4Run* aRun = nullptr;
  for ( auto& ru : *this ) {
      auto anotherRun = ru->GenerateRun();
      if ( aRun != nullptr && anotherRun != nullptr ) {
          G4Exception("G4MultiRunAction::GenerateRun()","Run0036",FatalException,
              "More than one registered UserRunAction return an instance"\
              " of G4Run, not allowed.");
          return nullptr;
      }
      aRun = anotherRun;
  }
  return aRun;
}

void G4MultiRunAction::BeginOfRunAction(const G4Run* run) {
  std::for_each( begin() , end() ,
      [run](G4UserRunActionUPtr& e) { e->BeginOfRunAction(run); }
  );
}

void G4MultiRunAction::EndOfRunAction(const G4Run* run) {
  std::for_each( begin() , end() ,
      [run](G4UserRunActionUPtr& e) { e->EndOfRunAction(run); }
  );
}

void G4MultiRunAction::SetMaster(G4bool val) {
  G4UserRunAction::SetMaster(val);
  std::for_each( begin() , end() ,
      [val](G4UserRunActionUPtr& e) { e->SetMaster(val); }
  );

}

