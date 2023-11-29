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

// Author: Ivana Hrivnacova, 10/09/2014  (ivana@ipno.in2p3.fr)

#include "G4BaseFileManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4Threading.hh"

//_____________________________________________________________________________
G4BaseFileManager::G4BaseFileManager(const G4AnalysisManagerState& state)
  : fState(state)
{}

//
// public methods
//

//_____________________________________________________________________________
void G4BaseFileManager::AddFileName(const G4String& fileName)
{
  // G4cout << "registering " << fileName << " in manager of " << GetFileType() << G4endl;

  // Do nothing in file name is already present
  for ( const auto& name : fFileNames ) {
    if ( name == fileName ) return;
  }

  fFileNames.push_back(fileName);
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetFileType() const
{
  return G4StrUtil::to_lower_copy(fState.GetType());
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetFullFileName(const G4String& baseFileName,
                                            G4bool isPerThread) const
{
  G4String fileName(baseFileName);
  if ( fileName == "" ) fileName = fFileName;

  // Take out file extension
  auto name = G4Analysis::GetBaseName(fileName);

  // Add thread Id to a file name if MT processing
  if ( isPerThread && ! fState.GetIsMaster() ) {
    std::ostringstream os;
    os << G4Threading::G4GetThreadId();
    name.append("_t");
    name.append(os.str());
  }

  // Add (back if it was present or is defined) file extension
  auto extension = G4Analysis::GetExtension(fileName, GetFileType());
  if (extension.size() != 0u) {
    name.append(".");
    name.append(extension);
  }

  return name;
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetHnFileName(const G4String& hnType,
                                          const G4String& hnName) const
{
  return G4Analysis::GetHnFileName(fFileName, GetFileType(), hnType, hnName);
}

//_____________________________________________________________________________
G4String  G4BaseFileManager::GetHnFileName(const G4String& fileName,
                                           G4int cycle) const
{
  // Do nothing if cycle is supported by the output type
  if (HasCycles()) return fileName;

  return G4Analysis::GetHnFileName(fileName, GetFileType(), cycle);
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetNtupleFileName(const G4String& ntupleName,
                                              G4int cycle) const
{
  // Do not pass cycle if supported by the output type
  auto cycleToPass = (HasCycles()) ? 0 : cycle;

  return G4Analysis::GetNtupleFileName(fFileName, GetFileType(), ntupleName, cycleToPass);
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetNtupleFileName(G4int ntupleFileNumber,
                                              G4int cycle) const
{
  // Do not pass cycle if supported by the output type
  auto cycleToPass = (HasCycles()) ? 0 : cycle;

  return G4Analysis::GetNtupleFileName(fFileName, GetFileType(), ntupleFileNumber, cycleToPass);
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetPlotFileName() const
{
  return G4Analysis::GetPlotFileName(fFileName);
}


