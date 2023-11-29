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

// /vis/geometry commands - John Allison  31st January 2006

#include "G4VisCommandsGeometry.hh"

#include "G4UIcmdWithAString.hh"
#include "G4VisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UImanager.hh"

#define G4warn G4cout

std::map<G4LogicalVolume*, const G4VisAttributes*>
G4VVisCommandGeometry::fVisAttsMap;

G4VVisCommandGeometry::~G4VVisCommandGeometry()
{
  // Delete all vis atts that were "new".  Do something like "restore"
  // without the "rebuild".
}

////////////// /vis/geometry/list ///////////////////////////////////////

G4VisCommandGeometryList::G4VisCommandGeometryList()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/geometry/list", this);
  fpCommand -> SetGuidance("Lists vis attributes of logical volume(s).");
  fpCommand -> SetGuidance("\"all\" lists all logical volumes.");
  fpCommand -> SetParameterName("logical-volume-name", omitable = true);
  fpCommand -> SetDefaultValue("all");
}

G4VisCommandGeometryList::~G4VisCommandGeometryList()
{
  delete fpCommand;
}

G4String G4VisCommandGeometryList::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometryList::SetNewValue(G4UIcommand*, G4String newValue)
{
  G4LogicalVolumeStore *pLVStore = G4LogicalVolumeStore::GetInstance();
  G4bool found = false;
  for (size_t iLV = 0; iLV < pLVStore->size(); iLV++ ) {
    G4LogicalVolume*pLV = (*pLVStore)[iLV];
    const G4String& logVolName = pLV->GetName();
    if (newValue == "all" || logVolName == newValue) {
      const G4VisAttributes* visAtts = pLV->GetVisAttributes();
      G4cout << "\nLogical Volume \"" << pLV->GetName() << "\":";
      if (visAtts) {
        G4cout << '\n' << *visAtts;
      } else {
        G4cout << " no vis attributes";
      }
      G4cout << G4endl;
    }
    if (logVolName == newValue) found = true;
  }
  if (newValue != "all" && !found) {
    if (fpVisManager->GetVerbosity() >= G4VisManager::errors) {
      G4warn << "ERROR: Logical volume \"" << newValue
	     << "\" not found in logical volume store." << G4endl;
    }
    return;
  }
}

////////////// /vis/geometry/restore ///////////////////////////////////////

G4VisCommandGeometryRestore::G4VisCommandGeometryRestore()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/geometry/restore", this);
  fpCommand -> SetGuidance("Restores vis attributes of logical volume(s).");
  fpCommand -> SetParameterName("logical-volume-name", omitable = true);
  fpCommand -> SetDefaultValue("all");
}

G4VisCommandGeometryRestore::~G4VisCommandGeometryRestore()
{
  delete fpCommand;
}

G4String G4VisCommandGeometryRestore::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometryRestore::SetNewValue(G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4LogicalVolumeStore *pLVStore = G4LogicalVolumeStore::GetInstance();
  size_t nLV = pLVStore->size();
  size_t iLV;
  G4LogicalVolume* pLV = 0;
  G4bool found = false;
  for (iLV = 0; iLV < nLV; iLV++ ) {
    pLV = (*pLVStore)[iLV];
    const G4String& logVolName = pLV->GetName();
    if (logVolName == newValue) found = true;
    if (newValue == "all" || logVolName == newValue) {
      VisAttsMapIterator i = fVisAttsMap.find(pLV);
      if (i != fVisAttsMap.end()) {
	const G4VisAttributes* newVisAtts = pLV->GetVisAttributes();
	const G4VisAttributes* oldVisAtts = i->second;
	pLV->SetVisAttributes(oldVisAtts);
	if (verbosity >= G4VisManager::confirmations) {
	  G4cout << "\nLogical Volume \"" << pLV->GetName()
		 << "\": re-setting vis attributes:\nwas: " << *newVisAtts
		 << "\nnow: " << *oldVisAtts
		 << G4endl;
	}
      }
    }
  }
  if (newValue != "all" && !found) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: Logical volume \"" << newValue
	     << "\" not found in logical volume store." << G4endl;
    }
    return;
  }
  if (fpVisManager->GetCurrentViewer()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/rebuild");
  }
}
