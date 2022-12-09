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

#include "G4XmlNtupleFileManager.hh"
#include "G4XmlFileManager.hh"
#include "G4XmlNtupleManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4AutoLock.hh"

using namespace G4Analysis;
using std::make_shared;

//_____________________________________________________________________________
G4XmlNtupleFileManager::G4XmlNtupleFileManager(const G4AnalysisManagerState& state)
 : G4VNtupleFileManager(state, "csv")
{}

//
// public methods
//

//_____________________________________________________________________________
std::shared_ptr<G4VNtupleManager> G4XmlNtupleFileManager::CreateNtupleManager()
{
  fNtupleManager = std::make_shared<G4XmlNtupleManager>(fState);
  fNtupleManager->SetFileManager(fFileManager);

  return fNtupleManager;
}

//_____________________________________________________________________________
G4bool G4XmlNtupleFileManager::ActionAtOpenFile(const G4String& /*fileName*/)
{
  // Create ntuples if they are booked
  // (The files will be created with creating ntuples)
  fNtupleManager->CreateNtuplesFromBooking(
    fBookingManager->GetNtupleBookingVector());

  return true;
}

//_____________________________________________________________________________
G4bool G4XmlNtupleFileManager::ActionAtWrite()
{
  auto ntupleVector = fNtupleManager->GetNtupleDescriptionVector();

  for ( auto ntupleDescription : ntupleVector ) {
    if (ntupleDescription->GetNtuple() != nullptr) {
      ntupleDescription->GetNtuple()->write_trailer();
    }
  }

  return true;
}

//_____________________________________________________________________________
G4bool G4XmlNtupleFileManager::ActionAtCloseFile()
{
  auto result = true;

  // Close ntuple files
  auto ntupleVector = fNtupleManager->GetNtupleDescriptionVector();
  for ( auto ntupleDescription : ntupleVector) {
    result &= fFileManager->CloseNtupleFile(ntupleDescription);
  }

  return result;
}

//_____________________________________________________________________________
G4bool G4XmlNtupleFileManager::Reset()
{
// Reset histograms and ntuple

  return fNtupleManager->Reset();
}

