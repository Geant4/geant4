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

#pragma once
#include "globals.hh"
#include "G4String.hh"
#include "G4MoleculeDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DamagedDeoxyribose : public G4MoleculeDefinition
{
    private:
        static G4DamagedDeoxyribose* fgInstance;
        G4DamagedDeoxyribose() {}
        ~G4DamagedDeoxyribose() {}

    public:
        static G4DamagedDeoxyribose* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DamagedAdenine : public G4MoleculeDefinition
{
    private:
        static G4DamagedAdenine* fgInstance;
        G4DamagedAdenine() {}
        ~G4DamagedAdenine() {}

    public:
        static G4DamagedAdenine* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DamagedGuanine : public G4MoleculeDefinition
{
    private:
        static G4DamagedGuanine* fgInstance;
        G4DamagedGuanine() {}
        ~G4DamagedGuanine() {}

    public:
        static G4DamagedGuanine* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DamagedThymine : public G4MoleculeDefinition
 {
 private:
     static G4DamagedThymine* fgInstance;
     G4DamagedThymine() {}
     ~G4DamagedThymine() {}

 public:
     static G4DamagedThymine* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DamagedCytosine : public G4MoleculeDefinition
{
    private:
        static G4DamagedCytosine* fgInstance;
        G4DamagedCytosine() {}
        ~G4DamagedCytosine() {}

    public:
        static G4DamagedCytosine* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Deoxyribose : public G4MoleculeDefinition
{
    private:
        static G4Deoxyribose* fgInstance;
        G4Deoxyribose() {}
        ~G4Deoxyribose() {}

    public:
        static G4Deoxyribose* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Phosphate : public G4MoleculeDefinition
{
    private:
        static G4Phosphate* fgInstance;
        G4Phosphate() {}
        ~G4Phosphate() {}

 public:
     static G4Phosphate* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Adenine : public G4MoleculeDefinition
 {
 private:
     static G4Adenine* fgInstance;
     G4Adenine() {}
     ~G4Adenine() {}

 public:
     static G4Adenine* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Guanine : public G4MoleculeDefinition
 {
 private:
     static G4Guanine* fgInstance;
     G4Guanine() {}
     ~G4Guanine() {}

 public:
     static G4Guanine* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Thymine : public G4MoleculeDefinition
 {
 private:
     static G4Thymine* fgInstance;
     G4Thymine() {}
     ~G4Thymine() {}

 public:
     static G4Thymine* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Cytosine : public G4MoleculeDefinition
{
private:
    static G4Cytosine* fgInstance;
    G4Cytosine() {}
    ~G4Cytosine() {}

public:
    static G4Cytosine* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4ModifiedHistone : public G4MoleculeDefinition
{
private:
    static G4ModifiedHistone* fgInstance;
    G4ModifiedHistone() {}
    ~G4ModifiedHistone() {}

public:
    static G4ModifiedHistone* Definition();
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Histone : public G4MoleculeDefinition
{
private:
    static G4Histone* fgInstance;
    G4Histone() {}
    ~G4Histone() {}

public:
    static G4Histone* Definition();
 };
