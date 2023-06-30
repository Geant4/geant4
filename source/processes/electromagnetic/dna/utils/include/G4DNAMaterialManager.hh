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
// 1/2/2023: Hoang: this file is used to check avaiable DNA materils,
// and keeps DNA cross sections.

#ifndef   G4DNAMaterialManager_hh
#define  G4DNAMaterialManager_hh 1

#include "globals.hh"
#include "G4DNACrossSectionDataSet.hh"
#include <map>
class G4Material;
class G4ParticleDefinition;
class G4StateManager;
class G4VEmModel;
enum class DNAModelType
{
  fDNAIonisation = 0,
  fDNAExcitation,
  fDNAElastics,
  fDNADefault
};

class G4DNAMaterialManager {
public:
  using MaterialMap = std::map<size_t, G4Material*>;

  static G4DNAMaterialManager *Instance();

  G4DNAMaterialManager(const G4DNAMaterialManager &) = delete;

  const G4DNAMaterialManager &operator=(const G4DNAMaterialManager &) = delete;

  G4VEmModel* GetModel(const DNAModelType& t);

  void SetMasterDataModel(const DNAModelType& t, G4VEmModel* m);

  G4bool IsLocked() const;

private:
  G4DNAMaterialManager();

  ~G4DNAMaterialManager() = default;

  static G4DNAMaterialManager *theInstance;

  MaterialMap fMaterials;
  G4StateManager*  fStateManager = nullptr;
  //this data stores only master models.
  //this data is only to read then non conflits
  std::map<DNAModelType,G4VEmModel*> fData;
};

#endif
