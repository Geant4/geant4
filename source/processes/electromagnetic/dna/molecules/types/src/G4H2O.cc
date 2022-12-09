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
// Author: Mathieu Karamitors 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4H2O.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

G4H2O* G4H2O::fgpInstance = nullptr;

G4H2O* G4H2O::Definition()
{
    if (fgpInstance != nullptr)
    {
        return fgpInstance;
    }

    const G4String name = "H2O";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if (pInstance == nullptr)
    {
        const G4String formatedName = "H_{2}O";

        G4double mass = 18.0153 * g / Avogadro * c_squared;
        pInstance = new G4MoleculeDefinition(name,
                                              mass,
                                              2e-5 * (cm * cm / s), // diffusion coefficient
                                              0, // charge
                                              8, // electronic levels
                                              2.75 * angstrom, // radius
                                              3, // atoms number
                                              0 // lifetime
        ); //picosecond set in dissociation process

        ((G4MoleculeDefinition*) pInstance)->SetLevelOccupation(0);
        ((G4MoleculeDefinition*) pInstance)->SetLevelOccupation(1);
        ((G4MoleculeDefinition*) pInstance)->SetLevelOccupation(2);
        ((G4MoleculeDefinition*) pInstance)->SetLevelOccupation(3);
        ((G4MoleculeDefinition*) pInstance)->SetLevelOccupation(4);
        ((G4MoleculeDefinition*) pInstance)->SetFormatedName(formatedName);

    }
    fgpInstance = static_cast<G4H2O*>(pInstance);
    return fgpInstance;
}

G4H2O* G4H2O::DefinitionIfExists()
{
    return fgpInstance;
}


