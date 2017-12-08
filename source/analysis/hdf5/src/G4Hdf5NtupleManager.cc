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
// $Id$

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#include "G4Hdf5NtupleManager.hh"
#include "G4Hdf5FileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4UnitsTable.hh"

#include "tools/ntuple_booking"

#include <iostream>
//#include <cstdio>

using namespace G4Analysis;

//_____________________________________________________________________________
G4Hdf5NtupleManager::G4Hdf5NtupleManager(const G4AnalysisManagerState& state)
 : G4TNtupleManager<tools::hdf5::ntuple>(state),
   fFileManager(nullptr)
{}

//_____________________________________________________________________________
G4Hdf5NtupleManager::~G4Hdf5NtupleManager()
{}

// 
// private methods
//
//_____________________________________________________________________________
void G4Hdf5NtupleManager::CreateTNtuple(
  G4TNtupleDescription<tools::hdf5::ntuple>* /*ntupleDescription*/,
  const G4String& /*name*/, const G4String& /*title*/)
{
  // Ntuple will be created at finish ntuple from ntuple_booking
}

//_____________________________________________________________________________
void G4Hdf5NtupleManager::CreateTNtupleFromBooking(
  G4TNtupleDescription<tools::hdf5::ntuple>* ntupleDescription)
{
    // create a file for this ntuple
    // if ( ! fFileManager->CreateNtupleFile(ntupleDescription) ) return;

    // Check ntuple directory
    if ( fFileManager->GetNtupleDirectory() < 0 ) {
      G4String inFunction = "G4Hdf5NtupleManager::::CreateTNtupleFromBooking";
      G4ExceptionDescription description;
      description << "      " 
        << "Cannot create ntuple. Ntuple directory does not exist." << G4endl;
      G4Exception(inFunction, "Analysis_W002", JustWarning, description);
      return;
    }

    auto basketSize = fFileManager->GetBasketSize();
    // auto compressionLevel = fState.GetCompressionLevel();
    auto compressionLevel = 0;

    // create ntuple
    ntupleDescription->fNtuple
      = new tools::hdf5::ntuple(
              G4cout, fFileManager->GetNtupleDirectory(), ntupleDescription->fNtupleBooking, 
              compressionLevel, basketSize);

    fNtupleVector.push_back(ntupleDescription->fNtuple);  
}

//_____________________________________________________________________________
void G4Hdf5NtupleManager::FinishTNtuple(
  G4TNtupleDescription<tools::hdf5::ntuple>* /*ntupleDescription*/)
{
  // if ( ! ntupleDescription->fNtuple ) return;

  // Only ntuple description is created on CreateNtuple call
  // CreateTNtupleFromBooking(ntupleDescription);

  fFileManager->LockNtupleDirectoryName();
}
