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
// G4tgrMaterialPropertiesTable
//
#include "G4tgrMaterialPropertiesTable.hh"
#include "G4tgrUtils.hh"

//---------------------------------------------------------------------
G4tgrMaterialPropertiesTable::G4tgrMaterialPropertiesTable(const std::vector<G4String>& wl)
{
  // :PROP name
  //   const_name1 value1
  //   array_name1 array_values1
  fwords = wl;
  //---------- Set name
  theName = G4tgrUtils::GetString(fwords[1]);
  G4cout << "Setting optical properties of " << theName << ":" << G4endl;
  G4tgrPopMaterialPropertiesTable(2);
}
  
//---------------------------------------------------------------------
G4tgrMaterialPropertiesTable::~G4tgrMaterialPropertiesTable()
{
}

//---------------------------------------------------------------------
void G4tgrMaterialPropertiesTable::G4tgrPopMaterialPropertiesTable(G4int idx)
{
    for (size_t i=idx; i<fwords.size(); i++) {
    G4String item_name = fwords[i];
    G4StrUtil::to_upper(item_name);
    if (G4StrUtil::contains(item_name, "SCINTILLATIONYIELD") 
      || item_name=="RESOLUTIONSCALE" 
      || G4StrUtil::contains(item_name, "TIMECONSTANT") 
      || G4StrUtil::contains(item_name, "MIEHG_") 
      || G4StrUtil::contains(item_name, "SCINTILLATIONRISETIME")) {
      theConstants.emplace_back(item_name, G4tgrUtils::GetDouble(fwords[i+1]));
      G4cout << "TextGeom: " << item_name << "=" << fwords[i+1] << G4endl;
      i++; // item_name value has been used
    } else if (item_name.substr(0,12)=="PHOTON_ENERG") {
      photonEnergyUnDefined=false;
      ArrayCount = G4tgrUtils::GetInt(fwords[i+1]);
      energies = new double[ArrayCount]; // create energy array
      
      for (int j=0; j<ArrayCount; j++)
        energies[j]=G4tgrUtils::GetDouble(fwords[i+2+j]);
      i=i+1+ArrayCount; // array has been used
    } else { // wavelength-dependent properties
      if (photonEnergyUnDefined) {
        G4cout<<"TextGeom: photon energies undefined, "
          <<"ignore all wavelength-dependent properties!"<<G4endl;
        break;
      }
      ArrayProperty ap; 
      ap.itemName = item_name;
      ap.arrayCount = ArrayCount;
      ap.values   = new double[ArrayCount];
      ap.energies = new double[ArrayCount];
      for (int j=0; j<ArrayCount; j++)
      {
        ap.values[j] = G4tgrUtils::GetDouble(fwords[i+1+j]);
        ap.energies[j] = energies[j];
      }
      theProperties.push_back(ap);
      i=i+ArrayCount;
    }
  }
}

  