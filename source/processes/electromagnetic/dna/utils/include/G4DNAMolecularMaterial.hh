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
// $Id: G4DNAMolecularMaterial.hh 103042 2017-03-10 11:50:07Z gcosmo $
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

/** 
 * \struct CompareMaterial
 * \brief Materials can be described as a derivation of existing "parent" 
 * materials in order to alter few of their features, such as density.
 * \p CompareMaterial compare materials taking into account
 * their possible "affiliation".
 */
struct CompareMaterial
{
  bool operator()(const G4Material* mat1, const G4Material* mat2) const;
};

typedef std::map<const G4Material*, double, CompareMaterial> ComponentMap;

/**
 * \class G4DNAMolecularMaterial
 * \brief G4DNAMolecularMaterial builds tables of molecular densities for chosen
 * molecular materials. The class handles homogeneous, composite and
 * derived materials. A material of interest is labeled as molecular if built
 * using the number of atoms rather than the mass fractions.
 *
 * \details
 * - Initialization:
 * G4DNAMolecularMaterial is initialized when 
 *    G4ApplicationState == G4State_Idle.
 * It should be initialized on the master thread and used in read-only mode 
 * during stepping. The singleton is thread-shared.
 *
 * - For Developers:
 * Use GetNumMolPerVolTableFor(molecule) in the concrete implementation of 
 *  G4VEmModel::Initialise or G4VProcess::PreparePhysicsTable
 * at run initialization to retrieve a read-only, thread-safe, table. 
 * The table is then built on the master thread at initialization time and 
 * shared between all threads and models.
 *
 * \note A G4material is labeled as molecular if built using the number of atoms
 *
 */

class G4DNAMolecularMaterial: public G4VStateDependent
{
public:
  static G4DNAMolecularMaterial* Instance();
  static void DeleteInstance();
  void Initialize();
  void Clear();

  virtual G4bool Notify(G4ApplicationState requestedState);
  
  //----------------------------------------------------------------------------

  /**
   * \fn const std::vector<double>* \
   *     GetDensityTableFor(const G4Material* searchedMaterial) const
   * \brief Retrieve a table of volumetric mass densities (mass per unit volume)
   * in the G4 unit system for chosen material.
   *
   * @param[in] searchedMaterial
   * The material which you'd like to retrieve the volumic mass
   * @pre The \p searchedMaterial used in parameter must be built as a
   * molecular material, using the number of atoms rather than the density
   * fractions.
   * \return
   * Pointer to a table of molecular densities for the \p searchedMaterial
   * indexed on the (parent) material index.
   *
   */
  const std::vector<double>* GetDensityTableFor(const G4Material*) const;
  
  /**
   * \fn const std::vector<double>* \
   *     GetNumMolPerVolTableFor(const G4Material* searchedMaterial) const
   * \brief Retrieve a table of molecular densities (number of molecules per
   * unit volume) in the G4 unit system for chosen material.
   *
   * @param[in] searchedMaterial
   * The material which you'd like to retrieve the molecular density
   * @pre The \p searchedMaterial used in parameter must be built as a 
   * molecular material, using the number of atoms rather than the density 
   * fractions.
   * \return
   * Pointer to a table of molecular densities for the \p searchedMaterial 
   * indexed on the (parent) material index.
   */
  const std::vector<double>* GetNumMolPerVolTableFor(const G4Material*) const;
  
  inline const std::vector<ComponentMap>* GetMassFractionTable() const{
    return fpCompFractionTable;
  }
  inline const std::vector<ComponentMap>* GetDensityTable() const{
    return fpCompDensityTable;
  }
  
  //----------------------------------------------------------------------------
  
  G4MolecularConfiguration* GetMolecularConfiguration(const G4Material*) const;
  
  /**
   * \fn void SetMolecularConfiguration(const G4Material* material, \
   *                                    G4MolecularConfiguration* molConf)
   * \brief Associate a molecular configuration to a G4material.
   *
   * @param[in] material
   * Pointer to a G4 material. The material
   * does not need to be defined as a molecular material.
   * @param[in] molConf
   * The molecular configuration corresponding to
   * the G4 \p material.
   */
  void SetMolecularConfiguration(const G4Material*,
                                 G4MolecularConfiguration*);
  
  /**
   * \fn void SetMolecularConfiguration(const G4Material* material, \
   *                                    const G4String& molConf)
   * \brief Associate a molecular configuration to a G4material.
   *
   * @param[in] material
   * Pointer to a G4 material. The material
   * does not need to be defined as a molecular material.
   * @param[in] molConf
   * User ID of the molecular configuration corresponding to
   * the G4 \p material.
   */
  void SetMolecularConfiguration(const G4Material*,
                                 const G4String&);
  
  /**
   * \fn void SetMolecularConfiguration(const G4Material* material, \
   *                                    const G4String& molConf)
   * \brief Associate a molecular configuration to a G4material.
   *
   * @param[in] material
   * Name of the G4 material. The material
   * does not need to be defined as a molecular material.
   * @param[in] molConf
   * User ID of the molecular configuration corresponding to
   * the G4 \p material.
   */
  void SetMolecularConfiguration(const G4String& materialName,
                                 const G4String& molUserIF);
  
  //----------------------------------------------------------------------------
  
  /**
   * \brief Deprecated
   * \deprecated Will return a G4 fatal exception.
   * Use instead GetNumMolPerVolTableFor(molecule) at run
   * initialization to retrieve a read-only, thread-safe, table.
   * \note A G4material is labeled as molecular if built using
   * the number of atoms.
   */
  G4double GetNumMoleculePerVolumeUnitForMaterial(const G4Material *mat);
  
  /**
   * \brief Deprecated
   * \deprecated Will return a G4 fatal exception.
   * Use instead GetNumMolPerVolTableFor(molecule) at run
   * initialization to retrieve a read-only, thread-safe, table.
   * \note A G4material is labeled as molecular if built using
   * the number of atoms.
   */
  G4double GetNumMolPerVolForComponentInComposite(const G4Material *composite,
                                                  const G4Material *component,
                                                  G4double massFraction);

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

  // Tables built for all molecular materials at initialization
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
