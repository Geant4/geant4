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
//
//

// G4SDStructure
#include "G4SDStructure.hh"
#include "G4ios.hh"

G4SDStructure::G4SDStructure(const G4String& aPath)
  : verboseLevel(0)
{
  pathName = aPath;
  dirName  = aPath;
  auto i  = dirName.length();
  if(i > 1)
  {
    dirName.erase(i - 1);
    auto isl = dirName.rfind('/');
    dirName.erase(0, isl + 1);
    dirName += "/";
  }
}

G4SDStructure::~G4SDStructure()
{
  for(auto st : structure)
    delete st;
  structure.clear();
  for(auto dt : detector)
    delete dt;
  detector.clear();
}

G4bool G4SDStructure::operator==(const G4SDStructure& right) const
{
  return (this == &right);
}

void G4SDStructure::AddNewDetector(G4VSensitiveDetector* aSD,
                                   const G4String& treeStructure)
{
  G4String remainingPath = treeStructure;
  remainingPath.erase(0, pathName.length());
  if(!remainingPath.empty())
  {  // The detector should be kept in subdirectory.
     // First, check if the subdirectoy exists.
    G4String subD         = ExtractDirName(remainingPath);
    G4SDStructure* tgtSDS = FindSubDirectory(subD);
    if(tgtSDS == nullptr)
    {  // Subdirectory not found. Create a new directory.
      subD.insert(0, pathName);
      tgtSDS = new G4SDStructure(subD);
      structure.push_back(tgtSDS);
    }
    tgtSDS->AddNewDetector(aSD, treeStructure);
  }
  else
  {  // The sensitive detector should be kept in this directory.
    G4VSensitiveDetector* tgtSD = GetSD(aSD->GetName());
    if(!tgtSD)
    {
      detector.push_back(aSD);
    }
    else if(tgtSD != aSD)
    {
#ifdef G4VERBOSE
      G4ExceptionDescription ed;
      ed << aSD->GetName() << " had already been stored in " << pathName
         << ". Object pointer is overwritten.\n";
      ed << "It's users' responsibility to delete the old sensitive detector "
            "object.";
      G4Exception("G4SDStructure::AddNewDetector()", "DET1010", JustWarning,
                  ed);
#endif
      RemoveSD(tgtSD);
      detector.push_back(aSD);
    }
  }
}

G4SDStructure* G4SDStructure::FindSubDirectory(const G4String& subD)
{
  for(auto st : structure)
  {
    if(subD == st->dirName)
      return st;
  }
  return nullptr;
}

G4VSensitiveDetector* G4SDStructure::GetSD(const G4String& aSDName)
{
  for(auto det : detector)
  {
    if(aSDName == det->GetName())
      return det;
  }
  return nullptr;
}

void G4SDStructure::RemoveSD(G4VSensitiveDetector* sd)
{
  auto det = std::find(detector.begin(), detector.end(), sd);
  if(det != detector.end())
    detector.erase(det);
}

G4String G4SDStructure::ExtractDirName(const G4String& aName)
{
  G4String subD = aName;
  auto i       = aName.find('/');
  if(i != G4String::npos)
    subD.erase(i + 1);
  return subD;
}

void G4SDStructure::Activate(const G4String& aName, G4bool sensitiveFlag)
{
  G4String aPath = aName;
  aPath.erase(0, pathName.length());
  if(aPath.find('/') != std::string::npos)
  {  // Command is ordered for a subdirectory.
    G4String subD         = ExtractDirName(aPath);
    G4SDStructure* tgtSDS = FindSubDirectory(subD);
    if(tgtSDS == nullptr)
    {  // The subdirectory is not found
      G4cout << subD << " is not found in " << pathName << G4endl;
    }
    else
    {
      tgtSDS->Activate(aName, sensitiveFlag);
    }
  }
  else if(aPath.empty())
  {  // Command is ordered for all detectors in this directory.
    for(auto det : detector)
      det->Activate(sensitiveFlag);
    for(auto st : structure)
      st->Activate(G4String("/"), sensitiveFlag);
  }
  else
  {  // Command is ordered to a particular detector.
    G4VSensitiveDetector* tgtSD = GetSD(aPath);
    if(tgtSD == nullptr)
    {  // The detector is not found.
      G4cout << aPath << " is not found in " << pathName << G4endl;
    }
    else
    {
      tgtSD->Activate(sensitiveFlag);
    }
  }
}

G4VSensitiveDetector* G4SDStructure::FindSensitiveDetector(
  const G4String& aName, G4bool warning)
{
  G4String aPath = aName;
  aPath.erase(0, pathName.length());
  if(aPath.find('/') != std::string::npos)
  {  // SD exists in sub-directory
    G4String subD         = ExtractDirName(aPath);
    G4SDStructure* tgtSDS = FindSubDirectory(subD);
    if(tgtSDS == nullptr)
    {  // The subdirectory is not found
      if(warning)
        G4cout << subD << " is not found in " << pathName << G4endl;
      return nullptr;
    }
    else
    {
      return tgtSDS->FindSensitiveDetector(aName, warning);
    }
  }
  else
  {  // SD must exist in this directory
    G4VSensitiveDetector* tgtSD = GetSD(aPath);
    if(tgtSD == nullptr)
    {  // The detector is not found.
      if(warning)
        G4cout << aPath << " is not found in " << pathName << G4endl;
    }
    return tgtSD;
  }
}

void G4SDStructure::Initialize(G4HCofThisEvent* HCE)
{
  // Broadcast to subdirectories.
  for(auto st : structure)
  {
    st->Initialize(HCE);
  }
  // Initialize all detectors in this directory.
  for(auto dt : detector)
  {
    if(dt->isActive())
      dt->Initialize(HCE);
  }
}

void G4SDStructure::Terminate(G4HCofThisEvent* HCE)
{
  // Broadcast to subdirectories.
  for(auto st : structure)
  {
    st->Terminate(HCE);
  }
  // Terminate all detectors in this directory.
  for(auto dt : detector)
  {
    if(dt->isActive())
      dt->EndOfEvent(HCE);
  }
}

void G4SDStructure::ListTree()
{
  G4cout << pathName << G4endl;
  for(auto sd : detector)
  {
    G4cout << pathName << sd->GetName();
    if(sd->isActive())
    {
      G4cout << "   *** Active ";
    }
    else
    {
      G4cout << "   XXX Inactive ";
    }
    G4cout << G4endl;
  }
  for(auto st : structure)
    st->ListTree();
}
