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

// Author: Ivana Hrivnacova, 15/09/2020  (ivana@ipno.in2p3.fr)

#include "G4VNtupleFileManager.hh"

#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include <utility>

using namespace G4Analysis;

namespace {

//_____________________________________________________________________________
void NtupleMergingWarning(std::string_view className,
                          std::string_view functionName,
                          const G4String& outputType)
{
  Warn("Ntuple merging is not available with " + outputType + " output.\n" +
       "Setting is ignored.", className, functionName);
}

}

//_____________________________________________________________________________
G4VNtupleFileManager::G4VNtupleFileManager(const G4AnalysisManagerState& state, G4String fileType)
  : fState(state),
    fFileType(std::move(fileType))
{}

///_____________________________________________________________________________
void G4VNtupleFileManager::SetNtupleMerging(G4bool /*mergeNtuples*/,
                   G4int /*nofReducedNtupleFiles*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetNtupleMerging", fFileType);
}

//_____________________________________________________________________________
void G4VNtupleFileManager::SetNtupleRowWise(G4bool /*rowWise*/,
                                          G4bool /*rowMode*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetNtupleRowWise", fFileType);
}

//_____________________________________________________________________________
void G4VNtupleFileManager::SetBasketSize(unsigned int /*basketSize*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetBasketSize", fFileType);
}

//_____________________________________________________________________________
void G4VNtupleFileManager::SetBasketEntries(unsigned int /*basketEntries*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetBasketEntries", fFileType);
}
