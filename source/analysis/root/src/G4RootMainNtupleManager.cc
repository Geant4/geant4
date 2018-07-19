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
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#include "G4RootMainNtupleManager.hh"
#include "G4RootNtupleManager.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/ntuple"

//_____________________________________________________________________________
G4RootMainNtupleManager::G4RootMainNtupleManager(G4RootNtupleManager* ntupleBuilder,
                                                 G4bool rowWise,
                                                 const G4AnalysisManagerState& state)
 : G4BaseAnalysisManager(state),
   fNtupleBuilder(ntupleBuilder),
   fRowWise(rowWise),
   fNtupleDirectory(nullptr),
   fNtupleVector()
{}

//_____________________________________________________________________________
G4RootMainNtupleManager::~G4RootMainNtupleManager()
{
  // ntuple objects are deleted automatically when closing a file 
}

//
// protected functions
//

//_____________________________________________________________________________
void G4RootMainNtupleManager::CreateNtuple(const tools::ntuple_booking& ntupleBooking,
                                           G4bool warn)
{
// Create ntuple from booking if file was open

  // Check that file is set
  if ( ! fNtupleDirectory ) {
    if ( warn ) {
      G4ExceptionDescription description;
      description 
        << "      " << "Ntuple file must be defined first." 
        << G4endl 
        << "      " << "Cannot create main ntuples from builder.";
        G4Exception("G4RootAnalysisManager::CreateNtuplesFromBooking",
                  "Analysis_W002", JustWarning, description);
      }
    return;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "main ntuple", ntupleBooking.name());
#endif

  // Create ntuple
  auto ntuple = new tools::wroot::ntuple(*fNtupleDirectory, ntupleBooking, fRowWise);
         // ntuple object is deleted automatically when closing a file
  auto basketSize = fNtupleBuilder->GetBasketSize();
  ntuple->set_basket_size(basketSize);

  fNtupleVector.push_back(ntuple);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()
      ->Message("create", "main ntuple", ntupleBooking.name());
#endif
}

//_____________________________________________________________________________
void G4RootMainNtupleManager::CreateNtuplesFromBooking()
{
// Create ntuple from booking in master 

  // Check that file is set
  if ( ! fNtupleDirectory ) {
    G4ExceptionDescription description;
    description 
      << "      " << "Ntuple file must be defined first." 
      << G4endl 
      << "      " << "Cannot create main ntuples from builder.";
      G4Exception("G4RootAnalysisManager::CreateNtuplesFromBooking",
                "Analysis_W002", JustWarning, description);
    return;
  }

  auto& ntupleDescriptionVector 
    = fNtupleBuilder->GetNtupleDescriptionVector();

  for ( auto ntupleDescription : ntupleDescriptionVector ) {    
    CreateNtuple(ntupleDescription->fNtupleBooking);
  }
}   

//_____________________________________________________________________________
G4bool G4RootMainNtupleManager::Merge()
{
  for ( auto ntuple : fNtupleVector ) {
    ntuple->merge_number_of_entries();
  }

  return true;
}

//_____________________________________________________________________________
G4bool G4RootMainNtupleManager::Reset(G4bool deleteNtuple)
{
  for ( auto ntuple : fNtupleVector ) {
    if ( deleteNtuple ) {
      delete ntuple;
    }  
  }

  fNtupleVector.clear(); 
  
  return true;
}
