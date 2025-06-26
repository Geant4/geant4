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

// Author: Christian Velten (2025)

#ifndef G4MOLECULELOCATOR_HH
#define G4MOLECULELOCATOR_HH 1
#pragma once

#include "G4TouchableHandle.hh"
#include "G4ThreadLocalSingleton.hh"

#include <memory>

class G4Track;
class G4ITNavigator;

class G4MoleculeLocator final
{
  friend class G4ThreadLocalSingleton<G4MoleculeLocator>;
  public:
    static G4MoleculeLocator* Instance();

    ~G4MoleculeLocator() = default;

  private:
    G4MoleculeLocator();
    G4MoleculeLocator(const G4MoleculeLocator&) = delete;
    G4MoleculeLocator(G4MoleculeLocator&&) = delete;
    G4MoleculeLocator& operator=(const G4MoleculeLocator&) = delete;
    G4MoleculeLocator& operator=(G4MoleculeLocator&&) = delete;

    G4ThreadLocalStatic G4MoleculeLocator* fpInstance;

    G4bool fIsInitialized{false};
    std::unique_ptr<G4ITNavigator> fNavigator;

    void Initialize();

  public:
    void LocateMoleculeSetStateAndTouchable(G4Track*);
    G4TouchableHandle LocateMoleculeTrack(const G4Track*);
};

#endif
