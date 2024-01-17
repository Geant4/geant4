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
/// \file UserMolecule.hh
/// \brief Definition of the UserMolecule class

#ifndef UserMolecule_h
#define UserMolecule_h

#include "G4Molecule.hh"
class UserMolecule : public G4Molecule
{
public:
    using G4Molecule::G4Molecule;
    ~UserMolecule() override = default;
    //From G4VUserTrackInformation
    void Print() const override {;;;};
    //  new/delete operators are overloded to use G4Allocator
    inline void *operator new(size_t);
#ifdef __IBMCPP__
    inline void *operator new(size_t sz, void* p)
    {
        return p;
    }
#endif
    inline void operator delete(void*);
      // Copy number
    void SetCopyNumber(G4int copyNum) {fCopyNumber=copyNum;}
    G4int GetCopyNumber() const {return fCopyNumber;}
    // DNA strand
    void SetStrand(G4int strand) {fStrand=strand;}
    G4int GetStrand() const {return fStrand;}

private:
    G4int fCopyNumber=-1;
    G4int fStrand=-1;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#if defined G4EM_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<UserMolecule>*& UserMoleculeAllocator();
#else
extern G4DLLIMPORT G4Allocator<UserMolecule>*& UserMoleculeAllocator();
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void * UserMolecule::operator new(size_t)
{
    if (!UserMoleculeAllocator())
    {
        UserMoleculeAllocator() = new G4Allocator<UserMolecule>;
    }
    return (void *)UserMoleculeAllocator()->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void UserMolecule::operator delete(void * aMolecule)
{
    UserMoleculeAllocator()->FreeSingle((UserMolecule *)aMolecule);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif