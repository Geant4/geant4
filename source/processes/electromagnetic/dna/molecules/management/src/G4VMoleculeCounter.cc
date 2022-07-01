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
//  G4VMoleculeCounter.cc
//  Geant4
//
//  Created by Mathieu Karamitros on 02/11/2016.
//
//

#include "G4VMoleculeCounter.hh"
#include "G4MoleculeCounter.hh"
G4ThreadLocal
G4VMoleculeCounter* G4VMoleculeCounter::fpInstance = nullptr;

//------------------------------------------------------------------------------

void
G4VMoleculeCounter::SetInstance(G4VMoleculeCounter* pCounterInstance)
{
    if (fpInstance != nullptr)
    {
        G4ExceptionDescription errMsg;
        errMsg << "The G4MoleculeCounter was already initialized." << G4endl
               << "The previous instance will be deleted in order to use yours." << G4endl
               << "However this can generate conflicts. Make sure you call G4MoleculeCounter::SetInstance"
                  "at the beginning of your application."
               << "A good place would be ActionInitialization::Build & BuildForMaster"
               << G4endl;

        G4Exception("G4MoleculeCounter::SetInstance",
                    "SINGLETON_ALREADY_INITIALIZED",
                    JustWarning, errMsg);
        delete fpInstance;
        fpInstance = nullptr;
    }

    fpInstance = pCounterInstance;
}

//------------------------------------------------------------------------------

G4VMoleculeCounter* G4VMoleculeCounter::Instance()
{
    if (fpInstance == nullptr)
    {
        fpInstance = new G4MoleculeCounter();
    }
    return fpInstance;
}

//------------------------------------------------------------------------------

void G4VMoleculeCounter::DeleteInstance()
{
    if (fpInstance != nullptr)
    {
        delete fpInstance;
        fpInstance = nullptr;
    }
}

//------------------------------------------------------------------------------

void G4VMoleculeCounter::InitializeInstance()
{
    if (fpInstance)
    {
        fpInstance->Initialize();
    }
}

//------------------------------------------------------------------------------

void G4VMoleculeCounter::Use(G4bool flag)
{
    fUse = flag;
}

//------------------------------------------------------------------------------

G4bool G4VMoleculeCounter::InUse()
{
    return fUse;
}

