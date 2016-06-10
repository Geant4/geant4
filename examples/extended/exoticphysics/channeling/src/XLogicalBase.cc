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

#include "XLogicalBase.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

XLogicalBase::XLogicalBase(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalBase::XLogicalBase(G4Element* vElement,XLogicalAtomicLattice* vLattice){
    SetElement(vElement);
    SetLattice(vLattice);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
XLogicalBase::~XLogicalBase(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalAtomicLattice* XLogicalBase::GetLattice(){
    return fLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Element* XLogicalBase::GetElement(){
    return fElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalBase::SetLattice(XLogicalAtomicLattice* vLattice){
    fLattice = vLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalBase::SetElement(G4Element* vElement){
    fElement = vElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XLogicalBase::ComputeAtomicFormFactor(){
    return 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4complex XLogicalBase::
ComputeStructureFactorSingleAtomicKind(G4int h,
                                       G4int k,
                                       G4int l){
    G4double vAtomicFormFactor = ComputeAtomicFormFactor();
    G4complex vResult =
        GetLattice()->ComputeGeometricalStructureFactorSingleKind(h,k,l);
    vResult =
        G4complex(vResult.real() * vAtomicFormFactor,
                  vResult.imag() * vAtomicFormFactor);
    return vResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


