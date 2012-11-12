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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4MOLECULEHANDLEMANAGER_HH
#define G4MOLECULEHANDLEMANAGER_HH

#include "globals.hh"
#include "G4ReferenceCountedHandle.hh"
#include <map>
#include <CLHEP/Utility/memory.h>

class G4Molecule;

typedef CLHEP::shared_ptr<const G4Molecule> G4MoleculeHandle;

class G4MoleculeHandleManager
{
public:
    static G4MoleculeHandleManager* Instance();
    ~G4MoleculeHandleManager();
    G4MoleculeHandle GetMoleculeHandle(const G4Molecule*);

    static void DeleteInstance();

private:
    G4MoleculeHandleManager();
    static G4MoleculeHandleManager* fInstance;

    struct CompMoleculePointer
    {
        G4bool operator()(const G4Molecule* mol1, const G4Molecule* mol2) const;
    };

    typedef std::map<const G4Molecule*, G4MoleculeHandle, G4MoleculeHandleManager::CompMoleculePointer> MoleculeHandleMap;
    MoleculeHandleMap fMoleculeHandle;
};

#endif // G4MOLECULEHANDLEMANAGER_HH
