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
// $Id: G4DNAChemistryManager.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4DNACHEMISTRYMANAGER_HH
#define G4DNACHEMISTRYMANAGER_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <memory>

class G4Track;
class G4DNAWaterExcitationStructure;
class G4DNAWaterIonisationStructure;
class G4Molecule;

enum ElectronicModification
{
    eIonizedMolecule,
    eExcitedMolecule,
    eDissociativeAttachment
};

/**
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

class G4DNAChemistryManager
{
    friend class std::auto_ptr<G4DNAChemistryManager>;
    ~G4DNAChemistryManager();

public:
    static G4DNAChemistryManager* Instance();

    /**
      * You should rather use DeleteInstance than the destructor of this class
      */
    static void DeleteInstance();

    /**
      * Tells the chemMan to write into a file
      * the position and electronic state of the water molecule
      * and the position thermalized or not of the solvated electron
      */
    void WriteInto(const G4String&, std::ios_base::openmode mode = std::ios_base::out);

    /** Close the file specified with WriteInto
      */
    void CloseFile();
    inline G4bool IsChemistryActived();
    inline void SetChemistryActivation(G4bool);

    /**
      * Method used by DNA physics model to create a water molecule.
      * The ElectronicModification is a flag telling wheter the molecule
      * is ionized or excited, the electronic level is calculated by the
      * model and the IncomingTrack is the track responsible for the creation
      * of this molecule, for instance an electron.
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
      * Before calling this method, it is also possible to check if the chemistry is activated
      * through IsChemistryActived().
      * This method will create the track corresponding to the transfered molecule and will be in charge
      * of loading the new track to the system.
      */

    void PushMolecule(G4Molecule*& molecule,
                      double time, const G4ThreeVector& position, int parentID);

    /**
      * WARNING : In case chemistry is not activated, PushMoleculeAtParentTimeAndPlace
      * will take care of deleting the transfered molecule.
      * Before calling this method, it is also possible to check if the chemistry is activated
      * through IsChemistryActived().
      * This method will create the track corresponding to the transfered molecule and will be in charge
      * of loading the new track to the system.
      */
    void PushMoleculeAtParentTimeAndPlace(G4Molecule*& molecule,
                                          const G4Track* /*theIncomingTrack*/);
protected :
    G4DNAWaterExcitationStructure* GetExcitationLevel();
    G4DNAWaterIonisationStructure* GetIonisationLevel();

private:
    G4DNAChemistryManager();
    static std::auto_ptr<G4DNAChemistryManager> fInstance;
    bool fActiveChemistry;

    std::ofstream  fOutput;
    G4bool fWriteFile;

    G4DNAWaterExcitationStructure* fExcitationLevel;
    G4DNAWaterIonisationStructure* fIonisationLevel;
};

inline G4bool G4DNAChemistryManager::IsChemistryActived()
{
    return fActiveChemistry;
}

inline void G4DNAChemistryManager::SetChemistryActivation(G4bool flag)
{
    fActiveChemistry = flag;
}

#endif // G4DNACHEMISTRYMANAGER_HH
