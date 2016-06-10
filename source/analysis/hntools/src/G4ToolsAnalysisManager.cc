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
// $Id: G4ToolsAnalysisManager.cc 91116 2015-06-20 12:33:45Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4ToolsAnalysisManager.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4H3ToolsManager.hh"
#include "G4P1ToolsManager.hh"
#include "G4P2ToolsManager.hh"
#include "G4PlotManager.hh"
#include "G4MPIToolsManager.hh"

//_____________________________________________________________________________
G4ToolsAnalysisManager::G4ToolsAnalysisManager(const G4String& type, G4bool isMaster)
 : G4VAnalysisManager(type, isMaster),
   fH1Manager(nullptr),
   fH2Manager(nullptr),
   fH3Manager(nullptr),
   fP1Manager(nullptr),
   fP2Manager(nullptr)
{
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
}

//_____________________________________________________________________________
G4ToolsAnalysisManager::~G4ToolsAnalysisManager()
{}

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::PlotImpl() 
{

  // Only master thread performs plotting
  if ( G4Threading::IsWorkerThread() )  return true;

  auto finalResult = true;

  // Create plotter
  G4PlotManager plotManager(fState);
  plotManager.OpenFile(fVFileManager->GetPlotFileName());

  // H1
  auto result 
    = plotManager.PlotAndWrite<tools::histo::h1d>(fH1Manager->GetH1Vector(), 
                                                  fH1Manager->GetHnVector());
  finalResult = finalResult && result;

  // H2
  result 
    = plotManager.PlotAndWrite<tools::histo::h2d>(fH2Manager->GetH2Vector(), 
                                                  fH2Manager->GetHnVector());
  finalResult = finalResult && result;

  // H3
  // not yet available in tools

  // P1
  result 
    = plotManager.PlotAndWrite<tools::histo::p1d>(fP1Manager->GetP1Vector(), 
                                                  fP1Manager->GetHnVector());
  finalResult = finalResult && result;

  // P2
  // not yet available in tools

  result = plotManager.CloseFile();
  finalResult = finalResult && result;

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::MergeImpl(tools::histo::hmpi* hmpi) 
{

  // if ( G4Threading::IsWorkerThread() )  return true;

  if ( ! hmpi )  return false;

  G4bool finalResult = true;

  // Create MPI manager
  G4MPIToolsManager mpiToolsManager(fState, hmpi);

  // H1
  G4bool result 
    = mpiToolsManager.Merge<tools::histo::h1d>(fH1Manager->GetH1Vector(), 
                                               fH1Manager->GetHnVector());
  finalResult = finalResult && result;

  // H2
  result 
    = mpiToolsManager.Merge<tools::histo::h2d>(fH2Manager->GetH2Vector(), 
                                               fH2Manager->GetHnVector());
  finalResult = finalResult && result;

  // H3
  result 
    = mpiToolsManager.Merge<tools::histo::h3d>(fH3Manager->GetH3Vector(), 
                                               fH3Manager->GetHnVector());
  finalResult = finalResult && result;

  // P1
  result 
    = mpiToolsManager.Merge<tools::histo::p1d>(fP1Manager->GetP1Vector(), 
                                               fP1Manager->GetHnVector());
  finalResult = finalResult && result;

  // P2
  result 
    = mpiToolsManager.Merge<tools::histo::p2d>(fP2Manager->GetP2Vector(), 
                                               fP2Manager->GetHnVector());
  finalResult = finalResult && result;

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisManager::Reset()
{
// Reset histograms and profiles

  auto finalResult = true;
  
  auto result = fH1Manager->Reset();
  finalResult = finalResult && result;

  result = fH2Manager->Reset();
  finalResult = finalResult && result;
  
  result = fH3Manager->Reset();
  finalResult = finalResult && result;
  
  result = fP1Manager->Reset();
  finalResult = finalResult && result;
  
  result = fP2Manager->Reset();
  finalResult = finalResult && result;
  
  return finalResult;
}  
 

