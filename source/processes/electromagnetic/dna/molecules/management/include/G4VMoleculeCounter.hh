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
//  G4VMoleculeCounter.hh
//  Geant4
//
//  Created by Mathieu Karamitros on 02/11/2016.
//
//

#pragma once

#include <G4Types.hh>
#include <G4ios.hh>
#include "G4ThreeVector.hh"

class G4MolecularConfiguration;

class G4MoleculeDefinition;

class G4VMoleculeCounter
{
protected:
    static G4ThreadLocal G4VMoleculeCounter* fpInstance;
    G4bool fUse = false;

    G4VMoleculeCounter() = default;

    virtual ~G4VMoleculeCounter() = default;

public:
    static void SetInstance(G4VMoleculeCounter*);

    static void DeleteInstance();

    using Reactant = const G4MolecularConfiguration;

    /*
     * If no instance of G4VMoleculeCounter is provided
     * then use G4MoleculeCounter
     */
    static G4VMoleculeCounter* Instance();

    static void InitializeInstance();

    /*
     * If the molecule counter is used, it will be called
     * at every creation/deletion of a molecule to
     * to increase/decrease the number at a given time.
     */
    void Use(G4bool flag = true);

    G4bool InUse();

    //----------------------------------------------------

    virtual void Initialize() = 0;

    virtual void ResetCounter() = 0;

    virtual void AddAMoleculeAtTime(Reactant*,
                                    G4double time,
                                    const G4ThreeVector* position = nullptr,
                                    int number = 1) = 0;

    virtual void RemoveAMoleculeAtTime(Reactant*,
                                       G4double time,
                                       const G4ThreeVector* position = nullptr,
                                       int number = 1) = 0;

    /* The dynamics of the given molecule won't be recorded. */
    virtual void DontRegister(const G4MoleculeDefinition*)
    {
    }

    virtual bool IsRegistered(const G4MoleculeDefinition*)
    {
        return false;
    }

    virtual void RegisterAll()
    {
    }
};