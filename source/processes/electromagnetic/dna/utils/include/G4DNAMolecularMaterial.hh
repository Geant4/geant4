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
// $Id: G4DNAMolecularMaterial.hh 101354 2016-11-15 08:27:51Z gcosmo $
//
// Author: Mathieu Karamitros
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4DNAMolecularMaterial_HH
#define G4DNAMolecularMaterial_HH

#include "globals.hh"
#include "G4ios.hh"
#include <map>
#include <vector>
#include "G4VStateDependent.hh"

class G4Material;
class G4MolecularConfiguration;

struct CompareMaterial
{
  // If the materials derives from a base material,
  // it should be able to find the derived material using the base material.
  bool operator()(const G4Material* mat1, const G4Material* mat2) const;
};

typedef std::map<const G4Material*, double, CompareMaterial> ComponentMap;

// G4DNAMolecularMaterial is initialized when G4ApplicationState == G4State_Idle

class G4DNAMolecularMaterial : public G4VStateDependent
{
public:
  static G4DNAMolecularMaterial* Instance();
  static void DeleteInstance();
  void Initialize();
  void Clear();

  virtual G4bool Notify(G4ApplicationState requestedState);
  
  //----------------------------------------------------------------------------

  const std::vector<double>* GetDensityTableFor(const G4Material*) const;
  const std::vector<double>* GetNumMolPerVolTableFor(const G4Material*) const;
  
  inline const std::vector<ComponentMap>* GetMassFractionTable() const{
    return fpCompFractionTable;
  }
  inline const std::vector<ComponentMap>* GetDensityTable() const{
    return fpCompDensityTable;
  }
  
  //----------------------------------------------------------------------------
  
  G4MolecularConfiguration* GetMolecularConfiguration(const G4Material*) const;
  void SetMolecularConfiguration(const G4Material*,
                                 G4MolecularConfiguration*);
  void SetMolecularConfiguration(const G4Material*,
                                 const G4String&);
  
  void SetMolecularConfiguration(const G4String& materialName,
                                 const G4String& molUserIF);

protected:
  static G4DNAMolecularMaterial* fInstance;
  G4DNAMolecularMaterial();
  G4DNAMolecularMaterial(const G4DNAMolecularMaterial& right);
  G4DNAMolecularMaterial& operator=(const G4DNAMolecularMaterial&);
  virtual ~G4DNAMolecularMaterial();
  void Create();
  void InitializeNumMolPerVol();
  void InitializeDensity();
  void RecordMolecularMaterial(G4Material* parentMaterial,
                               G4Material* molecularMaterial,
                               G4double fraction);
  void SearchMolecularMaterial(G4Material* parentMaterial,
                               G4Material* material,
                               double currentFraction);

  void AddMaterial(const G4Material*, double fraction);

  void PrintNotAMolecularMaterial(const char* methodName,
                                  const G4Material* lookForMaterial) const;

  std::vector<ComponentMap>* fpCompFractionTable;
  std::vector<ComponentMap>* fpCompDensityTable;
  std::vector<ComponentMap>* fpCompNumMolPerVolTable;

  mutable std::map<const G4Material*, std::vector<double>*, CompareMaterial>
            fAskedDensityTable;
  mutable std::map<const G4Material*, std::vector<double>*, CompareMaterial>
            fAskedNumPerVolTable;
  mutable std::map<const G4Material*, bool, CompareMaterial> fWarningPrinted;
  
  std::map<int /*Material ID*/,
           G4MolecularConfiguration*> fMaterialToMolecularConf;

  G4bool fIsInitialized;
  size_t fNMaterials;
};

#endif // G4DNAMolecularMaterial_HH
