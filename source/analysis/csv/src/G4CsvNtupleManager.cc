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
// $Id: G4CsvNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4CsvNtupleManager.hh"
#include "G4CsvFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "G4Threading.hh"

#include <iostream>

using namespace G4Analysis;

// 
// utility methods
//

//_____________________________________________________________________________
G4CsvNtupleManager::G4CsvNtupleManager(const G4AnalysisManagerState& state)
 : G4TNtupleManager<tools::wcsv::ntuple>(state),
   fFileManager(nullptr),
   fIsCommentedHeader(true),
   fIsHippoHeader(false)
{}

//_____________________________________________________________________________
G4CsvNtupleManager::~G4CsvNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
void G4CsvNtupleManager::CreateTNtuple(
  G4TNtupleDescription<tools::wcsv::ntuple>* ntupleDescription,
  const G4String& /*name*/, const G4String& title)
{
  // Create ntuple if the file is open (what means here that
  // a filename was already set)
  if ( fFileManager->GetFileName().size() ) {
    if ( fFileManager->CreateNtupleFile(ntupleDescription) ) {
      ntupleDescription->fNtuple 
        = new tools::wcsv::ntuple(*(ntupleDescription->fFile));
           // ntuple object is deleted when closing a file
      (ntupleDescription->fNtuple)->set_title(title); 
      fNtupleVector.push_back(ntupleDescription->fNtuple);       
    }       
  }  
}

//_____________________________________________________________________________
void G4CsvNtupleManager::CreateTNtupleFromBooking(
  G4TNtupleDescription<tools::wcsv::ntuple>* ntupleDescription)
{
    // create a file for this ntuple
    if ( ! fFileManager->CreateNtupleFile(ntupleDescription) ) return;

    // create ntuple
    ntupleDescription->fNtuple
      = new tools::wcsv::ntuple(
              *(ntupleDescription->fFile), G4cerr, ntupleDescription->fNtupleBooking);
    fNtupleVector.push_back(ntupleDescription->fNtuple);    
 }

//_____________________________________________________________________________
void G4CsvNtupleManager::FinishTNtuple(
  G4TNtupleDescription<tools::wcsv::ntuple>* ntupleDescription)
{
  if ( ! ntupleDescription->fNtuple ) return;

  // Write header if ntuple already exists
  if ( ! WriteHeader(ntupleDescription->fNtuple) ) {
     G4ExceptionDescription description;
     description << "      " 
                 << "Writing ntuple header has failed. ";
     G4Exception("G4CsvNtupleManager::FinishTNtuple()",
                 "Analysis_W022", JustWarning, description);
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
