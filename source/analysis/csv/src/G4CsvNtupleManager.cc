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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4CsvNtupleManager.hh"
#include "G4CsvFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

using namespace G4Analysis;

//
// utility methods
//

//_____________________________________________________________________________
G4CsvNtupleManager::G4CsvNtupleManager(const G4AnalysisManagerState& state)
 : G4TNtupleManager<tools::wcsv::ntuple, std::ofstream>(state)
{}

//
// private methods
//

//_____________________________________________________________________________
void G4CsvNtupleManager::CreateTNtupleFromBooking(
  CsvNtupleDescription* ntupleDescription)
{
  // create a file for this ntuple
  if ( ! fFileManager->CreateNtupleFile(ntupleDescription) ) return;

  // create ntuple
  ntupleDescription->SetNtuple(
    new tools::wcsv::ntuple(
          *(ntupleDescription->GetFile()), G4cerr, ntupleDescription->GetNtupleBooking()));
  fNtupleVector.push_back(ntupleDescription->GetNtuple());
 }

//_____________________________________________________________________________
void G4CsvNtupleManager::FinishTNtuple(
  CsvNtupleDescription* ntupleDescription,
  G4bool /*fromBooking*/)
{

  // Do nothing if the base file name was not yet defined
  if (fFileManager->GetFileName().size() == 0u) return;

  // Create ntuple from booking
  if (ntupleDescription->GetNtuple() == nullptr) {
    CreateTNtupleFromBooking(ntupleDescription);
  }

  // Return if creating ntuple failed
  if (ntupleDescription->GetNtuple() == nullptr) {
    Warn("Creating ntuple has failed.", fkClass, "FinishTNtuple");
    return;
  }

  // Write header if ntuple already exists
  if ( ! WriteHeader(ntupleDescription->GetNtuple()) ) {
    Warn("Writing ntuple header has failed.", fkClass, "FinishTNtuple");
  }
}

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::WriteHeader(tools::wcsv::ntuple* ntuple) const
{
// Write header if ntuple already exists and if this option is activated.
// When both Hippo and Commented headers are selected, only Commented
// header, which reading is supported.
// Return false only if an error occurred.

  if ( fIsCommentedHeader ) {
    return ntuple->write_commented_header(G4cout);
  }

  // write hippo header (if activated and if not commented header)
  if ( fIsHippoHeader ) {
    ntuple->write_hippo_header();
    return true;
  }

  return true;
}
