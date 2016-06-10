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
// $Id: G4ToolsAnalysisReader.cc 91116 2015-06-20 12:33:45Z ihrivnac $

// Author: Ivana Hrivnacova, 20/07/2015  (ivana@ipno.in2p3.fr)

#include "G4ToolsAnalysisReader.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4H3ToolsManager.hh"
#include "G4P1ToolsManager.hh"
#include "G4P2ToolsManager.hh"
#include "G4MPIToolsManager.hh"

//_____________________________________________________________________________
G4ToolsAnalysisReader::G4ToolsAnalysisReader(const G4String& type, G4bool isMaster)
 : G4VAnalysisReader(type, isMaster),
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
G4ToolsAnalysisReader::~G4ToolsAnalysisReader()
{}

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4ToolsAnalysisReader::Reset()
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
