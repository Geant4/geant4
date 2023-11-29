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

// Author: Ivana Hrivnacova, 20/07/2015  (ivana@ipno.in2p3.fr)

#include "G4ToolsAnalysisReader.hh"

//_____________________________________________________________________________
G4ToolsAnalysisReader::G4ToolsAnalysisReader(const G4String& type)
 : G4VAnalysisReader(type)
{
  // Create managers
  fH1Manager = new G4THnToolsManager<1, tools::histo::h1d>(fState);
  fH2Manager = new G4THnToolsManager<2, tools::histo::h2d>(fState);
  fH3Manager = new G4THnToolsManager<3, tools::histo::h3d>(fState);
  fP1Manager = new G4THnToolsManager<2, tools::histo::p1d>(fState);
  fP2Manager = new G4THnToolsManager<3, tools::histo::p2d>(fState);
      // The managers will be deleted by the base class

  // Set managers to base class which takes then their ownership
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
  SetH3Manager(fH3Manager);
  SetP1Manager(fP1Manager);
  SetP2Manager(fP2Manager);
}

//
// protected methods
//

//_____________________________________________________________________________
G4int G4ToolsAnalysisReader::ReadH1Impl(const G4String& h1Name,
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool isUserFileName)
{
  return ReadTImpl<tools::histo::h1d>(
           h1Name, fileName, dirName, isUserFileName, fH1Manager);
}

//_____________________________________________________________________________
G4int G4ToolsAnalysisReader::ReadH2Impl(const G4String& h2Name,
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool isUserFileName)
{
  return ReadTImpl<tools::histo::h2d>(
           h2Name, fileName, dirName, isUserFileName, fH2Manager);
}

//_____________________________________________________________________________
G4int G4ToolsAnalysisReader::ReadH3Impl(const G4String& h3Name,
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool isUserFileName)
{
  return ReadTImpl<tools::histo::h3d>(
           h3Name, fileName, dirName, isUserFileName, fH3Manager);
}

//_____________________________________________________________________________
G4int G4ToolsAnalysisReader::ReadP1Impl(const G4String& p1Name,
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool isUserFileName)
{
  return ReadTImpl<tools::histo::p1d>(
           p1Name, fileName, dirName, isUserFileName, fP1Manager);
}

//_____________________________________________________________________________
G4int G4ToolsAnalysisReader::ReadP2Impl(const G4String& p2Name,
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool isUserFileName)
{
  return ReadTImpl<tools::histo::p2d>(
           p2Name, fileName, dirName, isUserFileName, fP2Manager);
}

//_____________________________________________________________________________
G4bool G4ToolsAnalysisReader::Reset()
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
