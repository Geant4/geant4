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
// $Id: G4RootRNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#include "G4RootRNtupleManager.hh"
#include "G4AnalysisManagerState.hh"

//_____________________________________________________________________________
G4RootRNtupleManager::G4RootRNtupleManager(const G4AnalysisManagerState& state)
 : G4TRNtupleManager<tools::rroot::ntuple>(state)
{}

//_____________________________________________________________________________
G4RootRNtupleManager::~G4RootRNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::GetTNtupleRow(
  G4TRNtupleDescription<tools::rroot::ntuple>* ntupleDescription)
{
  auto ntuple = ntupleDescription->fNtuple;
  auto ntupleBinding = ntupleDescription->fNtupleBinding;

  G4bool isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {

    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Ntuple initialization failed !!"; 
      G4Exception("G4RootRNtuple::GetTNtupleRow()",
                  "Analysis_WR021", JustWarning, description);
      return false;
    }
    ntupleDescription->fIsInitialized = true;
    ntuple->start();
  }

  auto next = ntuple->next();
  if ( next ) {
    if ( ! ntuple->get_row() ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Ntuple get_row() failed !!"; 
      G4Exception("G4RootRNtuple::GetTNtupleRow()",
                  "Analysis_WR021", JustWarning, description);
      return false;
    }
  }  

  return next;
}   
