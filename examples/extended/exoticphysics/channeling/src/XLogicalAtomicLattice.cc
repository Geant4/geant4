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

#include "XLogicalAtomicLattice.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

XLogicalAtomicLattice::XLogicalAtomicLattice(){
    InitializeXLogicalAtomicLattice();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalAtomicLattice::~XLogicalAtomicLattice(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLattice::InitializeXLogicalAtomicLattice(){
    fLatticeAtomNumber = 1;
    for(G4int i=0;i<MAXLATTICEATOMS;i++)
        fLatticeAtomPosition[i] = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XLogicalAtomicLattice::GetAtomPosition(G4int i){
    if(i<fLatticeAtomNumber){
        return fLatticeAtomPosition[i];
    }
    else{
        G4cout << "XLogicalAtomicLattice::GetAtomPosition - atom " <<
        i << " does not exist!!" <<std::endl;
    }
    return G4ThreeVector(-1.,-1.,-1.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int XLogicalAtomicLattice::GetLatticeNumberOfAtoms(){
    return fLatticeAtomNumber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLattice::AddAtom(G4ThreeVector vAtomPosition){
    fLatticeAtomNumber++;
    //Add an atom to the lattice
    fLatticeAtomPosition[fLatticeAtomNumber - 1] = vAtomPosition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalAtomicLattice::DeleteAtom(G4ThreeVector vAtomPosition){
    //Delete atoms in the lattice in the selected position
    
    for(G4int i=0;i<fLatticeAtomNumber;i++)
        if(vAtomPosition == fLatticeAtomPosition[i])
        {
            for(G4int j=(i+1);j<fLatticeAtomNumber;j++)
            {
                fLatticeAtomPosition[j-1]=fLatticeAtomPosition[j];
            }
            i--;
            fLatticeAtomNumber--;
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4complex XLogicalAtomicLattice::
ComputeGeometricalStructureFactorSingleKind(G4int h,
                                            G4int k,
                                            G4int l){
    G4double vTempDouble = 0.;
    G4complex vResult = G4complex(0.,0.);

    for(G4int i=0;i<fLatticeAtomNumber;i++)
    {
        vTempDouble = 0.0;
        vTempDouble += h * fLatticeAtomPosition[i].x();
        vTempDouble += k * fLatticeAtomPosition[i].y();
        vTempDouble += l * fLatticeAtomPosition[i].z();
        vResult += G4complex(std::cos(2 * CLHEP::pi * vTempDouble),
                             std::sin(2 * CLHEP::pi * vTempDouble));
    }

    return vResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
