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

#include "G4ToolsAnalysisManager.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4H3ToolsManager.hh"
#include "G4P1ToolsManager.hh"
#include "G4P2ToolsManager.hh"
#include "G4PlotManager.hh"
#include "G4MPIToolsManager.hh"
#include "G4AnalysisUtilities.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"

using namespace G4Analysis;

namespace {
  //Mutex to lock master manager when merging histograms
  G4Mutex mergeHnMutex = G4MUTEX_INITIALIZER;
}

//_____________________________________________________________________________
G4ToolsAnalysisManager* G4ToolsAnalysisManager::Instance()
{
  return fgToolsInstance;
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::IsInstance()
{
  return ( fgToolsInstance != nullptr );
}

//_____________________________________________________________________________
G4ToolsAnalysisManager::G4ToolsAnalysisManager(const G4String& type)
 : G4VAnalysisManager(type)
{
  // Set instance pointer
  if ( ! G4Threading::IsWorkerThread() ) fgMasterToolsInstance = this;
  fgToolsInstance = this;

  // Create managers
  fH1Manager = new G4H1ToolsManager(fState);
  fH2Manager = new G4H2ToolsManager(fState);
  fH3Manager = new G4H3ToolsManager(fState);
  fP1Manager = new G4P1ToolsManager(fState);
  fP2Manager = new G4P2ToolsManager(fState);
      // The managers will be deleted by the base class

  // Set managers to base class which takes then their ownership
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
  SetH3Manager(fH3Manager);
  SetP1Manager(fP1Manager);
  SetP2Manager(fP2Manager);

  // Plot manager
  SetPlotManager(std::make_shared<G4PlotManager>(fState));

  // Messenger
  fMessenger = std::make_unique<G4ToolsAnalysisMessenger>(this);
}

//_____________________________________________________________________________
G4ToolsAnalysisManager::~G4ToolsAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterToolsInstance = nullptr;
  fgToolsInstance = nullptr;
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::PlotImpl()
{

  // Only master thread performs plotting
  if ( G4Threading::IsWorkerThread() )  return true;

  auto result = true;

  // Open output file
  fPlotManager->OpenFile(fVFileManager->GetPlotFileName());

  // H1
  result
    &= fPlotManager->PlotAndWrite<tools::histo::h1d>(fH1Manager->GetH1Vector(),
                                                    fH1Manager->GetHnVector());

  // H2
  result
    &= fPlotManager->PlotAndWrite<tools::histo::h2d>(fH2Manager->GetH2Vector(),
                                                    fH2Manager->GetHnVector());

  // H3
  // not yet available in tools

  // P1
  result
    &= fPlotManager->PlotAndWrite<tools::histo::p1d>(fP1Manager->GetP1Vector(),
                                                    fP1Manager->GetHnVector());

  // P2
  // not yet available in tools

  // Close output file
  result &= fPlotManager->CloseFile();

  return result;
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::MergeImpl(tools::histo::hmpi* hmpi)
{

  // if ( G4Threading::IsWorkerThread() )  return true;

  if ( ! hmpi )  return false;

  auto result = true;

  // Create MPI manager
  G4MPIToolsManager mpiToolsManager(fState, hmpi);

  // H1
  result
    &= mpiToolsManager.Merge<tools::histo::h1d>(fH1Manager->GetH1Vector(),
                                               fH1Manager->GetHnVector());

  // H2
  result
    &= mpiToolsManager.Merge<tools::histo::h2d>(fH2Manager->GetH2Vector(),
                                               fH2Manager->GetHnVector());

  // H3
  result
    &= mpiToolsManager.Merge<tools::histo::h3d>(fH3Manager->GetH3Vector(),
                                               fH3Manager->GetHnVector());

  // P1
  result
    &= mpiToolsManager.Merge<tools::histo::p1d>(fP1Manager->GetP1Vector(),
                                               fP1Manager->GetHnVector());

  // P2
  result
    &= mpiToolsManager.Merge<tools::histo::p2d>(fP2Manager->GetP2Vector(),
                                               fP2Manager->GetHnVector());

  return result;
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::WriteImpl()
{
  // Nothing to be done on worker
  if ( G4Threading::IsWorkerThread() ) return false;

  auto result = true;

  // Write all histograms/profile on master
  result &=WriteT(fH1Manager->GetH1Vector(), fH1Manager->GetHnVector());
  result &=WriteT(fH2Manager->GetH2Vector(), fH2Manager->GetHnVector());
  result &=WriteT(fH3Manager->GetH3Vector(), fH3Manager->GetHnVector());
  result &=WriteT(fP1Manager->GetP1Vector(), fP1Manager->GetHnVector());
  result &=WriteT(fP2Manager->GetP2Vector(), fP2Manager->GetHnVector());

  return result;
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::ResetImpl()
{
// Reset histograms and profiles

  auto result = true;

  result &= fH1Manager->Reset();
  result &= fH2Manager->Reset();
  result &= fH3Manager->Reset();
  result &= fP1Manager->Reset();
  result &= fP2Manager->Reset();

  return result;
}

//_____________________________________________________________________________
void G4ToolsAnalysisManager::ClearImpl()
{
// Reset histograms and profiles

  fH1Manager->ClearData();
  fH2Manager->ClearData();
  fH3Manager->ClearData();
  fP1Manager->ClearData();
  fP2Manager->ClearData();
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::Merge()
{
  // Nothing to be done on master
  if ( ! G4Threading::IsWorkerThread() ) return false;

  if ( ! fgMasterToolsInstance ) {
    if (! IsEmpty() ) {
      Warn("No master G4AnalysisManager instance exists.\n"
           "Histogram/profile data will not be merged.",
           fkClass, "Merge");
      return false;
    }
    return true;
  }

  Message(kVL4, "merge on worker", "histograms");

  // The worker manager just adds its histograms to the master
  fH1Manager->Merge(mergeHnMutex, fgMasterToolsInstance->fH1Manager);
  fH2Manager->Merge(mergeHnMutex, fgMasterToolsInstance->fH2Manager);
  fH3Manager->Merge(mergeHnMutex, fgMasterToolsInstance->fH3Manager);
  fP1Manager->Merge(mergeHnMutex, fgMasterToolsInstance->fP1Manager);
  fP2Manager->Merge(mergeHnMutex, fgMasterToolsInstance->fP2Manager);

  Message(kVL3, "merge on worker", "histograms");

  return true;
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::IsEmpty()
{
  return fH1Manager->IsEmpty() && fH2Manager->IsEmpty() && fH3Manager->IsEmpty() &&
         fP1Manager->IsEmpty() && fP2Manager->IsEmpty();
}
