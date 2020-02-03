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

#pragma once

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
class G4UIcmdWithAnInteger;
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
    ~G4DNAChemistryManager() override;

public:
    //============================================================================
    // STATIC METHODS
    //============================================================================
    static G4DNAChemistryManager* Instance();
    static G4DNAChemistryManager* GetInstanceIfExists();

    static void DeleteInstance();
    static G4bool IsActivated();

    //============================================================================
    // VIRTUAL METHODS
    //============================================================================
    // G4VStateDependent
    G4bool Notify(G4ApplicationState requestedState) override;
    // G4UImessenger
    void SetNewValue(G4UIcommand*, G4String) override;
    G4String GetCurrentValue(G4UIcommand* pCommand) override;

    //============================================================================
    // INITIALIZATION AND FINALIZATION METHODS
    //============================================================================
    G4bool IsChemistryActivated();
    void SetChemistryActivation(G4bool);

    /** Chemistry list is managed outside the chemistry manager (eg. constructor). */
    void SetChemistryList(G4VUserChemistryList&);

    /** Not a constructor or when used in standalone? Prefer this method. */
    void SetChemistryList(std::unique_ptr<G4VUserChemistryList>);

    // [[deprecated]] : chemistry list should never be nullptr
    void SetChemistryList(G4VUserChemistryList*);

    void Deregister(G4VUserChemistryList&);

    void Initialize();
    void Run();
    void Clear();

    /**
     * @brief Inject custom species to the simulation
     * @details This method should be called per thread, possibly from
     * ActionInitialisation::Build.
     * One can decide to set the same gun for all threads.
     * It is the user responsibility to handle the pointer deletion.
     */
    void SetGun(G4ITGun* pChemSpeciesGun);

    void SetPhysChemIO(std::unique_ptr<G4VPhysChemIO> pPhysChemIO);

    void SetVerbose(G4int verbose);

    /**
     * If the chemistry module is used in standalone (ie. without running the physics
     * stage beforehand), the physics table still needs to be built.
     * It is therefore necessary to flag the chemistry module as being run
     * in standalone.
     */
    void UseAsStandalone(G4bool flag);
    G4bool IsCounterResetWhenRunEnds() const;

    void ResetCounterWhenRunEnds(G4bool resetCounterWhenRunEnds);

    void ForceMasterReinitialization();
    void TagThreadForReinitialization();
    void ForceThreadReinitialization();
    void ForceRebuildingPhysicsTable();

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
    void AddEmptyLineInOutputFile();

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
                             const G4Track* /*pIncomingTrack*/);

    /**
     * This method should be used by the physics model of the ElectronSolvatation
     * process.
     */
    void CreateSolvatedElectron(const G4Track* /*pIncomingTrack*/,
                                G4ThreeVector* pFinalPosition = nullptr);

    void PushMolecule(std::unique_ptr<G4Molecule> pMolecule,
                      G4double time,
                      const G4ThreeVector& position,
                      G4int parentID);

protected:
    void HandleStandaloneInitialization();
    void PushTrack(G4Track*);
    void SetGlobalTemperature(G4double temperatureKelvin);

    G4DNAWaterExcitationStructure* GetExcitationLevel();
    G4DNAWaterIonisationStructure* GetIonisationLevel();
    void InitializeFile();
    void InitializeMaster();
    void InitializeThread();
    void InitializeThreadSharedData();

    G4DNAChemistryManager();

private:
    std::unique_ptr<G4UIdirectory> fpChemDNADirectory;
    std::unique_ptr<G4UIcmdWithABool> fpActivateChem;
    std::unique_ptr<G4UIcmdWithAnInteger> fpRunChem;
    std::unique_ptr<G4UIcmdWithoutParameter> fpSkipReactionsFromChemList;
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpScaleForNewTemperature;
    std::unique_ptr<G4UIcmdWithoutParameter> fpInitChem;

    static G4DNAChemistryManager* fgInstance;
    G4bool fActiveChemistry;

    struct ThreadLocalData{
        ThreadLocalData();
        ~ThreadLocalData();
        std::unique_ptr<G4VPhysChemIO> fpPhysChemIO;
        G4bool fThreadInitialized = false;
    };

    static G4ThreadLocal ThreadLocalData* fpThreadData;

    G4bool fMasterInitialized;
    G4bool fForceThreadReinitialization;

    std::unique_ptr<G4DNAWaterExcitationStructure> fpExcitationLevel;
    std::unique_ptr<G4DNAWaterIonisationStructure> fpIonisationLevel;

    std::unique_ptr<G4VUserChemistryList> fpUserChemistryList;
    G4bool fOwnChemistryList;
    G4bool fUseInStandalone;
    G4bool fPhysicsTableBuilt;
    G4bool fSkipReactions;

    G4bool fGeometryClosed;

    G4int fVerbose;
    G4bool fResetCounterWhenRunEnds;
};
