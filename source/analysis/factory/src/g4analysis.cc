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

// Author: Ivana Hrivnacova, 09/09/2019  (ivana@ipno.in2p3.fr)

#include "g4analysis.hh"
#include "G4AnalysisUtilities.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4RootAnalysisManager.hh"
#include "G4XmlAnalysisManager.hh"
#ifdef TOOLS_USE_HDF5
#include "G4Hdf5AnalysisManager.hh"
#endif
#include "G4Exception.hh"

namespace {

void DoFatalException(const G4String& outputName)
{
  G4ExceptionDescription description;
  description
    << "    \"" << outputName << "\" output type is not supported." << G4endl
    << "    " << "Analysis manager cannot be created.";
  G4Exception("G4Analysis::ManagerInstance",
              "Analysis_F002", FatalException, description);  
}

}

namespace G4Analysis
{

G4ToolsAnalysisManager* ManagerInstance(const G4String& outputName)
{
  // Get output type
  G4AnalysisOutput outputType = GetOutput(outputName, false);
  if ( outputType == G4AnalysisOutput::kNone ) {
    DoFatalException(outputName);
    return nullptr;
  }

  // Create the analysis manager of a selected type
  switch ( outputType ) {
    case G4AnalysisOutput::kCsv:
      return G4CsvAnalysisManager::Instance();
      break;
#ifdef TOOLS_USE_HDF5
    case G4AnalysisOutput::kHdf5: 
      return G4Hdf5AnalysisManager::Instance();
      break;
#endif
    case G4AnalysisOutput::kRoot: 
      return G4RootAnalysisManager::Instance();
      break;
    case G4AnalysisOutput::kXml:
      return G4XmlAnalysisManager::Instance();
      break;
    case G4AnalysisOutput::kNone:
    default:
      break;
  }
  
  // should never reach this line
  DoFatalException(outputName);
  return nullptr;
}

}
