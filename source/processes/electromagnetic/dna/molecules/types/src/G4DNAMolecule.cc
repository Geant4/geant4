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

#include "G4DNAMolecule.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4MoleculeDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"

G4DamagedDeoxyribose* G4DamagedDeoxyribose::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DamagedDeoxyribose* G4DamagedDeoxyribose::Definition()
{
    const G4String name = "Damaged_Deoxyribose";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 134.1305 * g / Avogadro * c_squared;//wikipedia
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                               5, 3.0 * angstrom, // radius
                                               2 // number of atoms
                                               ); 
    }
    
    fgInstance = static_cast<G4DamagedDeoxyribose*>(pInstance);
    return fgInstance;
}

G4DamagedAdenine* G4DamagedAdenine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DamagedAdenine* G4DamagedAdenine::Definition()
{
     const G4String name = "Damaged_Adenine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
         G4double mass = 135.1267 * g / Avogadro * c_squared;
         pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                               5, 3.0 * angstrom, // radius
                                               2 // number of atoms
                                               );
    }
    fgInstance = static_cast<G4DamagedAdenine*>(pInstance);
    return fgInstance;
}


G4DamagedGuanine* G4DamagedGuanine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DamagedGuanine* G4DamagedGuanine::Definition()
{
    const G4String name = "Damaged_Guanine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
         G4double mass = 151.1261 * g / Avogadro * c_squared;
         pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                               5, 3.0 * angstrom, // radius
                                               2 // number of atoms
                                               );
    }
    fgInstance = static_cast<G4DamagedGuanine*>(pInstance);
    return fgInstance;
}

G4DamagedThymine* G4DamagedThymine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DamagedThymine* G4DamagedThymine::Definition()
{
    const G4String name = "Damaged_Thymine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
         G4double mass = 126.1133 * g / Avogadro * c_squared;
         pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                               5, 3.0 * angstrom, // radius
	                                           2 // number of atoms
                                               );
    }

    fgInstance = static_cast<G4DamagedThymine*>(pInstance);
    return fgInstance;
}

G4DamagedCytosine* G4DamagedCytosine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DamagedCytosine* G4DamagedCytosine::Definition()
{
    const G4String name = "Damaged_Cytosine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 111.102 * g / Avogadro * c_squared;
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                              5, 2.9 * angstrom, // radius
                                              2 // number of atoms
                                              );
     }
     fgInstance = static_cast<G4DamagedCytosine*>(pInstance);
     return fgInstance;
}


G4Deoxyribose* G4Deoxyribose::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Deoxyribose* G4Deoxyribose::Definition()
{
    const G4String name = "Deoxyribose";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
         G4double mass = 134.1305 * g / Avogadro * c_squared;
         pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                               5, 2.9 * angstrom, // radius
	                                           2 // number of atoms
                                               );
    }
    fgInstance = static_cast<G4Deoxyribose*>(pInstance);
    return fgInstance;
}

G4Phosphate* G4Phosphate::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Phosphate* G4Phosphate::Definition()
{
    const G4String name = "Phosphate";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 94.9714 * g / Avogadro * c_squared;
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                             5, 2.7 * angstrom, // radius
                                             2 // number of atoms
                                             );
    }
    fgInstance = static_cast<G4Phosphate*>(pInstance);
    return fgInstance;
}

G4Adenine* G4Adenine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Adenine* G4Adenine::Definition()
{
    const G4String name = "Adenine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 135.1267 * g / Avogadro * c_squared;
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                             5, 3 * angstrom, // radius
                                             1 // number of atoms
                                             );
    }
    fgInstance = static_cast<G4Adenine*>(pInstance);
    return fgInstance;
}

G4Guanine* G4Guanine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Guanine* G4Guanine::Definition()
{
    const G4String name = "Guanine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 151.1261 * g / Avogadro * c_squared;
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0, 
                                             5, 3 * angstrom, // radius
                                             1 // number of atoms
                                             );

    }

    fgInstance = static_cast<G4Guanine*>(pInstance);
    return fgInstance;
}

G4Thymine* G4Thymine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Thymine* G4Thymine::Definition()
{
    const G4String name = "Thymine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
    G4double mass = 126.1133 * g / Avogadro * c_squared;
    pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                         5, 3 * angstrom, // radius
                                         1 // number of atoms
                                         );
    }
    fgInstance = static_cast<G4Thymine*>(pInstance);
    return fgInstance;
}

G4Cytosine* G4Cytosine::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Cytosine* G4Cytosine::Definition()
{
    const G4String name = "Cytosine";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 111.102 * g / Avogadro * c_squared;
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                             5, 3 * angstrom, // radius
                                             1 // number of atoms
                                             );
    }
    fgInstance = static_cast<G4Cytosine*>(pInstance);
    return fgInstance;
}

G4ModifiedHistone* G4ModifiedHistone::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ModifiedHistone* G4ModifiedHistone::Definition()
{
    const G4String name = "Modified_Histone";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 1.4e4 * g / Avogadro * c_squared;
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                             5, 24 * angstrom, // radius
                                             1 // number of atoms
                                             );
    }
    fgInstance = static_cast<G4ModifiedHistone*>(pInstance);
    return fgInstance;
}


G4Histone* G4Histone::fgInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Histone* G4Histone::Definition()
{
    const G4String name = "Histone";
    if(fgInstance != nullptr)
    {
        return fgInstance;
    }
    G4ParticleTable* pTable =G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pInstance = pTable->FindParticle(name);

    if(pInstance == nullptr)
    {
        G4double mass = 1.4e4 * g / Avogadro * c_squared;//14 kDa in wikipedia
        pInstance = new G4MoleculeDefinition(name, mass, 0 * (m * m / s), 0,
                                             5, 24 * angstrom, // radius
                                             1 // number of atoms
                                             );
    }
    fgInstance = static_cast<G4Histone*>(pInstance);
    return fgInstance;
}
