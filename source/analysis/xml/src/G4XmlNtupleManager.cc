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
// $Id: G4XmlNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4XmlNtupleManager.hh"
#include "G4XmlFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4Threading.hh"
#include "G4UnitsTable.hh"

#include "tools/ntuple_booking"

#include <iostream>
//#include <cstdio>

using namespace G4Analysis;

//_____________________________________________________________________________
G4XmlNtupleManager::G4XmlNtupleManager(const G4AnalysisManagerState& state)
 : G4TNtupleManager<tools::waxml::ntuple>(state),
   fFileManager(nullptr)
{}

//_____________________________________________________________________________
G4XmlNtupleManager::~G4XmlNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
void G4XmlNtupleManager::CreateTNtuple(
  G4TNtupleDescription<tools::waxml::ntuple>* ntupleDescription,
  const G4String& /*name*/, const G4String& /*title*/)
{
  // Create ntuple if the file is open (what means here that
  // a filename was already set)
  if ( fFileManager->GetFileName().size() ) {
    if ( fFileManager->CreateNtupleFile(ntupleDescription) ) {
      ntupleDescription->fNtuple 
        = new tools::waxml::ntuple(*(ntupleDescription->fFile));
           // ntuple object is deleted when closing a file
      fNtupleVector.push_back(ntupleDescription->fNtuple);       
    }       
  }
}

//_____________________________________________________________________________
void G4XmlNtupleManager::CreateTNtupleFromBooking(
  G4TNtupleDescription<tools::waxml::ntuple>* ntupleDescription)
{
    // create a file for this ntuple
    if ( ! fFileManager->CreateNtupleFile(ntupleDescription) ) return;

    // create ntuple
    ntupleDescription->fNtuple
      = new tools::waxml::ntuple(
              *(ntupleDescription->fFile), G4cerr, ntupleDescription->fNtupleBooking);
    fNtupleVector.push_back(ntupleDescription->fNtuple);  
}

//_____________________________________________________________________________
void G4XmlNtupleManager::FinishTNtuple(
  G4TNtupleDescription<tools::waxml::ntuple>* ntupleDescription)
{
  if ( ! ntupleDescription->fNtuple ) return;

  G4String path = "/";
  path.append(fFileManager->GetNtupleDirectoryName());
  ntupleDescription->fNtuple
    ->write_header(path, ntupleDescription->fNtupleBooking.name(), 
                   ntupleDescription->fNtupleBooking.title());  
  fFileManager->LockNtupleDirectoryName();
}
