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

#include "XCrystalCharacteristicArray.hh"

XCrystalCharacteristicArray::XCrystalCharacteristicArray(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalCharacteristicArray::~XCrystalCharacteristicArray(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XCrystalCharacteristicArray::
ComputeEC(G4ThreeVector vPosition,
          XPhysicalLattice* vLattice){
    
    G4ThreeVector vValue = G4ThreeVector(0.,0.,0.);
    
    if(fCharacteristicVector.size()!=0){
        for(unsigned int i=0;i<fCharacteristicVector.size();i++){
            vValue += fCharacteristicVector.at(i)->GetEC(vPosition,vLattice);
        }
        return (vValue * G4double(1./G4double(fCharacteristicVector.size())));
    }
    
    return G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XCrystalCharacteristicArray::
ComputePositionInUnitCell(G4ThreeVector vPosition,
                          XPhysicalLattice*vLattice){
    
    G4ThreeVector vPositionInTheCellPrevious = G4ThreeVector(0.,0.,0.);
    G4ThreeVector vPositionInTheCell = G4ThreeVector(0.,0.,0.);
    
    if(fCharacteristicVector.size()!=0){
        for(unsigned int i=0;i<fCharacteristicVector.size();i++){
            if(vPositionInTheCellPrevious != vPositionInTheCell){
                return G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
            }
            vPositionInTheCellPrevious = vPositionInTheCell;
            vPositionInTheCell =
            fCharacteristicVector.at(i)->ComputePositionInUnitCell(vPosition,
                                                                   vLattice);
        }
        return vPositionInTheCell;
    }
    
    return G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

