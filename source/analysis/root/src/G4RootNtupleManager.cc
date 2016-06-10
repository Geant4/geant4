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
// $Id: G4RootNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4RootNtupleManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"

using namespace G4Analysis;

//_____________________________________________________________________________
G4RootNtupleManager::G4RootNtupleManager(const G4AnalysisManagerState& state)
 : G4TNtupleManager<tools::wroot::ntuple>(state),
   fNtupleDirectory(nullptr)
{}

//_____________________________________________________________________________
G4RootNtupleManager::~G4RootNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
void G4RootNtupleManager::CreateTNtuple(
  G4TNtupleDescription<tools::wroot::ntuple>* ntupleDescription,
  const G4String& name, const G4String& title)
{
  // Create ntuple if the file is open
  if ( fNtupleDirectory ) {
    ntupleDescription->fNtuple
      = new tools::wroot::ntuple(*fNtupleDirectory, name, title);
    ntupleDescription->fIsNtupleOwner = false;  
           // ntuple object is deleted automatically when closing a file
    fNtupleVector.push_back(ntupleDescription->fNtuple);       
  }
}

//_____________________________________________________________________________
void G4RootNtupleManager::CreateTNtupleFromBooking(
  G4TNtupleDescription<tools::wroot::ntuple>* ntupleDescription)
{
    ntupleDescription->fNtuple
      = new tools::wroot::ntuple(
              *fNtupleDirectory, ntupleDescription->fNtupleBooking);
    ntupleDescription->fIsNtupleOwner = false;  
           // ntuple object is deleted automatically when closing a file
    fNtupleVector.push_back(ntupleDescription->fNtuple);  
}

//_____________________________________________________________________________
void G4RootNtupleManager::FinishTNtuple(
  G4TNtupleDescription<tools::wroot::ntuple>* /*ntupleDescription*/)
{
  // nothing to be done here
}
