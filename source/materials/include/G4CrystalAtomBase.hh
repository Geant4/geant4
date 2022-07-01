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
//

//---------------------------------------------------------------------------
//
// ClassName:   G4CrystalAtomBase
//
// Description: Contains crystal properties
//
// Class description:
//
// XXX
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 21-04-16, created by E.Bagli

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4CrystalAtomBase_HH
#define G4CrystalAtomBase_HH 1

#include <utility>
#include <vector>
#include <G4ThreeVector.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4CrystalAtomBase
{
public:  // with description
    //
    // Constructor to create a lattice
    //
    G4CrystalAtomBase() {;};
    G4CrystalAtomBase(const G4ThreeVector& apos) { AddPos(apos); };
    ~G4CrystalAtomBase() {;};
    
    //
    // Atom positions in the lattice are stored in fractional coordinates
    // i.e. the position in the unit cell described as a fractional
    // position along each cell edge.
    // Atoms may be removed or added
    //
private:
    std::vector<G4ThreeVector> thePos;

public:
    inline
    std::vector<G4ThreeVector> GetPos() {return thePos;}

    inline
    G4ThreeVector GetPos(G4int idx) {return thePos[idx];}
    inline void AddPos(const G4ThreeVector& a3vec) { thePos.push_back(a3vec); }
    inline
    void SetPos(std::vector<G4ThreeVector> a3vecvec) {
      thePos = std::move(a3vecvec);
    }
    inline
    void DelPos(G4int idx) {thePos.erase(thePos.begin()+idx);}
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
