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

#include "G4Hdf5RNtupleManager.hh"

//_____________________________________________________________________________
G4Hdf5RNtupleManager::G4Hdf5RNtupleManager(const G4AnalysisManagerState& state)
 : G4TRNtupleManager<tools::hdf5::ntuple>(state)
{}

//_____________________________________________________________________________
G4Hdf5RNtupleManager::~G4Hdf5RNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4Hdf5RNtupleManager::GetTNtupleRow(
  G4TRNtupleDescription<tools::hdf5::ntuple>* ntupleDescription)
{
  auto ntuple = ntupleDescription->fNtuple;

  auto isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {
    auto ntupleBinding = ntupleDescription->fNtupleBinding;
    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Ntuple initialization failed !!"; 
      G4Exception("G4Hdf5RNtuple::GetNtupleRow()",
                  "Analysis_WR021", JustWarning, description);
      return false;
    }
    ntupleDescription->fIsInitialized = true;
  }

  if ( ! ntuple->get_row() ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Ntuple get_row() failed !!"; 
    G4Exception("G4Hdf5RNtuple::GetNtupleRow()",
                "Analysis_WR021", JustWarning, description);
    return false;
  }

  return true;
}   
