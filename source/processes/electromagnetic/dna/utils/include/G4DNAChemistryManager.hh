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
// $Id: G4DNAChemistryManager.hh 103042 2017-03-10 11:50:07Z gcosmo $
//

// Author: Mathieu Karamitros
//
// The code is developed in the framework of the ESA AO7146
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

#ifndef G4DNACHEMISTRYMANAGER_HH
#define G4DNACHEMISTRYMANAGER_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <memory>
#include "G4UImessenger.hh"
#include "G4VStateDependent.hh"

class G4Track;
class G4DNAWaterExcitationStructure;
class G4DNAWaterIonisationStructure;
class G4Molecule;
class G4VUserChemistryList;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4ITGun;
class G4VPhysChemIO;

enum ElectronicModification
{
  eIonizedMolecule,
  eExcitedMolecule,
  eDissociativeAttachment
};

/**
 * WARNING: THIS CLASS IS A PROTOTYPE
 * G4DNAChemistryManager is called from the physics models.
 * It creates the water molecules and the solvated electrons and
 * and send them to G4ITStepManager to be treated in the chemistry stage.
 * For this, the fActiveChemistry flag needs to be on.
 * It is also possible to give already molecule's pointers already built.
 * G4DNAChemistryManager will then be in charge of creating the track and loading
 * it to the IT system.
 * The user can also ask to create a file containing a information about the
 * creation of water molecules and solvated electrons.
 */

class G4DNAChemistryManager: public G4UImessenger,
                             public G4VStateDependent
{
protected:
  virtual ~G4DNAChemistryManager();

public:
  //============================================================================
  // STATIC METHODS
  //============================================================================
  static G4DNAChemistryManager* Instance();
  static G4DNAChemistryManager* GetInstanceIfExists();

  /**
   * You should rather use DeleteInstance than the destructor of this class
   */
  static void DeleteInstance();
  G4bool IsChemistryActivated();
  void SetChemistryActivation(G4bool);
  static G4bool IsActivated();
  static void Activated(G4bool flag = true);

  //============================================================================
  // VIRTUAL METHODS
  //============================================================================
  virtual G4bool Notify(G4ApplicationState requestedState);
  virtual void SetNewValue(G4UIcommand*, G4String);
  virtual G4String GetCurrentValue(G4UIcommand * command);

  //============================================================================
  // INITIALIZATION AND FINALIZATION METHODS
  //============================================================================
  void Initialize();
  inline void SetChemistryList(G4VUserChemistryList*);
  inline void Deregister(G4VUserChemistryList*);
  void SetGlobalTemperature(G4double temp_K);

  inline void ForceMasterReinitialization();
  inline void TagThreadForReinitialization();

  void Run();
  void Clear();
  void Gun(G4ITGun*, bool physicsTableToBuild = true);
  inline void ForceThreadReinitialization();
  inline void ForceRebuildingPhysicsTable();

  //============================================================================
  // FILE OPERATIONS
  //============================================================================
  /**
   * Tells the chemMan to write into a file
   * the position and electronic state of the water molecule
   * and the position thermalized or not of the solvated electron
   */
  void WriteInto(const G4String&, std::ios_base::openmode mode =
      std::ios_base::out);
  void AddEmptyLineInOuputFile();

  /**
   * Close the file specified with WriteInto
   */
  void CloseFile();

  //============================================================================
  // PUSH MOLECULES
  //============================================================================
  /**
   * Method used by DNA physics model to create a water molecule.
   * The ElectronicModification is a flag telling whether the molecule
   * is ionized or excited, the electronic level is calculated by the
   * model and the IncomingTrack is the track responsible for the creation
   * of this molecule (electron, proton...).
   */
  void CreateWaterMolecule(ElectronicModification,
                           G4int /*electronicLevel*/,
                           const G4Track* /*theIncomingTrack*/);

  /**
   * On the same idea as the previous method but for solvated electron.
   * This method should be used by the physics model of the ElectronSolvatation
   * process.
   */
  void CreateSolvatedElectron(const G4Track* /*theIncomingTrack*/,
                              G4ThreeVector* finalPosition = 0);

  /**
   * WARNING : In case chemistry is not activated, PushMolecule will take care
   * of deleting the transfered molecule.
   * Before calling this method, it is also possible to check if the chemistry
   * is activated through IsChemistryActived().
   * This method will create the track corresponding to the transfered molecule
   * and will be in charge of loading the new track to the system.
   */

  void PushMolecule(G4Molecule*& molecule,
                    G4double time,
                    const G4ThreeVector& position,
                    G4int parentID);

  /**
   * WARNING : In case chemistry is not activated, PushMoleculeAtParentTimeAndPlace
   * will take care of deleting the transfered molecule.
   * Before calling this method, it is also possible to check if the chemistry
   * is activated through IsChemistryActived().
   * This method will create the track corresponding to the transfered molecule
   * and will be in charge of loading the new track to the system.
   */
  void PushMoleculeAtParentTimeAndPlace(G4Molecule*& molecule,
                                        const G4Track* /*theIncomingTrack*/);


  inline void SetVerbose(G4int verbose)
  {
    fVerbose = verbose;
  }

  inline void SetBuildPhysicsTable(G4bool flag)
  {fBuildPhysicsTable = flag;}

  G4bool IsCounterResetWhenRunEnds() const
  {
    return fResetCounterWhenRunEnds;
  }

  void ResetCounterWhenRunEnds(G4bool resetCounterWhenRunEnds)
  {
    fResetCounterWhenRunEnds = resetCounterWhenRunEnds;
  }

  void SetPhysChemIO(G4VPhysChemIO* physChemIO);
  
protected:
  G4DNAWaterExcitationStructure* GetExcitationLevel();
  G4DNAWaterIonisationStructure* GetIonisationLevel();
  void InitializeFile();
  void InitializeMaster();
  void InitializeThread();

  G4DNAChemistryManager();

private:
  G4UIdirectory* fpChemDNADirectory;
  G4UIcmdWithABool* fpActivateChem;
  G4UIcmdWithoutParameter* fpRunChem;
  G4UIcmdWithoutParameter* fpSkipReactionsFromChemList;
  //G4UIcmdWithADoubleAndUnit* fpGridSize; // not used in release
  G4UIcmdWithADoubleAndUnit* fpScaleForNewTemperature;
  G4UIcmdWithoutParameter* fpInitChem;
  
  static G4DNAChemistryManager* fgInstance;
  G4bool fActiveChemistry;
  
  struct ThreadLocalData{
    ThreadLocalData();
    ~ThreadLocalData();
    G4VPhysChemIO* fpPhysChemIO;
    G4bool fThreadInitialized_tl;
  };
  
  static G4ThreadLocal ThreadLocalData* fpThreadData;
  
  G4bool fMasterInitialized;
  G4bool fForceThreadReinitialization;

  G4DNAWaterExcitationStructure* fpExcitationLevel;
  G4DNAWaterIonisationStructure* fpIonisationLevel;

  G4VUserChemistryList* fpUserChemistryList;
  G4bool fBuildPhysicsTable;
  G4bool fPhysicsTableBuilt;
  G4bool fSkipReactions;

  G4bool fGeometryClosed;

  G4int fVerbose;
  G4bool fResetCounterWhenRunEnds;
};

//------------------------------------------------------------------------------

inline void G4DNAChemistryManager::ForceRebuildingPhysicsTable()
{
  fPhysicsTableBuilt = false;
}

inline void G4DNAChemistryManager::SetChemistryList(G4VUserChemistryList* chemistryList)
{
  fpUserChemistryList = chemistryList;
  Activated();
}

inline void G4DNAChemistryManager::Deregister(G4VUserChemistryList* chemistryList)
{
  if (fpUserChemistryList == chemistryList) fpUserChemistryList = 0;
}

inline void G4DNAChemistryManager::ForceMasterReinitialization()
{
  fMasterInitialized = false;
  InitializeMaster();
}

inline void G4DNAChemistryManager::ForceThreadReinitialization()
{
  fForceThreadReinitialization = true;
}

inline void G4DNAChemistryManager::TagThreadForReinitialization()
{
  fpThreadData->fThreadInitialized_tl = false;
}

#endif // G4DNACHEMISTRYMANAGER_HH
