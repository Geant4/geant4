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

#include "XCrystalIntegratedDensityHub.hh"

XCrystalIntegratedDensityHub::XCrystalIntegratedDensityHub(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensityHub::~XCrystalIntegratedDensityHub(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::SetPotential(
                    XVCrystalCharacteristic* vPotential){
    fPotential = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XCrystalIntegratedDensityHub::GetPotential(){
    return fPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::SetDensityElectron(
                    XVCrystalCharacteristic* vDensity){
    fDensityElectron = vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XCrystalIntegratedDensityHub::GetDensityElectron(){
    return fDensityElectron;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::SetDensityNuclei(
                    XVCrystalCharacteristic* vDensity){
    fDensityNuclei = vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XCrystalIntegratedDensityHub::GetDensityNuclei(){
    return fDensityNuclei;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::SetXPhysicalLattice(
                    XPhysicalLattice* vLattice){
    fLattice = vLattice;
    fDensityElectron->InitializePhysicalLattice(fLattice);
    fDensityNuclei->InitializePhysicalLattice(fLattice);
    fPotential->InitializePhysicalLattice(fLattice);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice*  XCrystalIntegratedDensityHub::GetXPhysicalLattice(){
    return fLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::SetIntegratedDensityNuclei(
                                    XVCrystalIntegratedDensity* vDensity,
                                    G4int vParticleCharge){
    if(vParticleCharge < 0){
        fIntDensNucleiNegative = vDensity;
    }
    else if(vParticleCharge > 0){
            fIntDensNucleiPositive = vDensity;
    }
    else{
        G4cout << "XCrystalIntegratedDensityHub:: ERROR  - Charge == 0";
        G4cout << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalIntegratedDensity* XCrystalIntegratedDensityHub::
GetIntegratedDensityNuclei(G4int vParticleCharge){
    if(vParticleCharge < 0.){
        return fIntDensNucleiNegative;
    }
    else if(vParticleCharge > 0.){
        return fIntDensNucleiPositive;
    }
    else{
        G4cout << "XCrystalIntegratedDensityHub:: ERROR  - Charge == 0";
        G4cout << G4endl;
        return NULL;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::SetIntegratedDensityElectron(
                            XVCrystalIntegratedDensity* vDensity,
                            G4int vParticleCharge){
    if(vParticleCharge < 0){
        fIntDensElectronNegative = vDensity;
    }
    else if(vParticleCharge > 0){
        fIntDensElectronPositive = vDensity;
    }
    else{
        G4cout << "XCrystalIntegratedDensityHub:: ERROR  - Charge == 0";
        G4cout << G4endl;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalIntegratedDensity* XCrystalIntegratedDensityHub::
GetIntegratedDensityElectron(G4int vParticleCharge){
    if(vParticleCharge < 0.){
        return fIntDensElectronNegative;
    }
    else if(vParticleCharge > 0.){
        return fIntDensElectronPositive;
    }
    else{
        G4cout << "XCrystalIntegratedDensityHub:: ERROR  - Charge == 0" ;
        G4cout << G4endl;
        return NULL;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XCrystalIntegratedDensityHub::HasBeenInitialized(
                                    XPhysicalLattice* vLattice){
    //now it checks only of the table is initialized,
    //it does not check if the particular crystal is initialized.
    //To be changed in the future!
    if(vLattice != GetXPhysicalLattice()) return false;
    else return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::InitializeTables(){
    fIntDensElectronNegative->SetPotential(fPotential);
    fIntDensElectronNegative->SetDensity(fDensityElectron);
    fIntDensElectronNegative->SetXPhysicalLattice(fLattice);
    fIntDensElectronNegative->SetParticleCharge(-1);
    fIntDensElectronNegative->InitializeTable();

    fIntDensElectronPositive->SetPotential(fPotential);
    fIntDensElectronPositive->SetDensity(fDensityElectron);
    fIntDensElectronPositive->SetXPhysicalLattice(fLattice);
    fIntDensElectronPositive->SetParticleCharge(+1);
    fIntDensElectronPositive->InitializeTable();
    
    fIntDensNucleiNegative->SetPotential(fPotential);
    fIntDensNucleiNegative->SetDensity(fDensityNuclei);
    fIntDensNucleiNegative->SetXPhysicalLattice(fLattice);
    fIntDensNucleiNegative->SetParticleCharge(-1);
    fIntDensNucleiNegative->InitializeTable();
    
    fIntDensNucleiPositive->SetPotential(fPotential);
    fIntDensNucleiPositive->SetDensity(fDensityNuclei);
    fIntDensNucleiPositive->SetXPhysicalLattice(fLattice);
    fIntDensNucleiPositive->SetParticleCharge(+1);
    fIntDensNucleiPositive->InitializeTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalIntegratedDensityHub::
GetIntegratedDensityNuclei(G4double vPotential,
                           XPhysicalLattice* vLattice,
                           G4int vParticleCharge){
    if(vParticleCharge < 0.){
        return fIntDensNucleiNegative->GetIntegratedDensity(vPotential,
                                                            vLattice,
                                                            vParticleCharge);
    }
    else if(vParticleCharge > 0.){
        return fIntDensNucleiPositive->GetIntegratedDensity(vPotential,
                                                            vLattice,
                                                            vParticleCharge);
    }
    else{
        G4cout << "XCrystalIntegratedDensityHub:: ERROR  - Charge == 0";
        G4cout << G4endl;
        return -1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalIntegratedDensityHub::
GetIntegratedDensityElectron(G4double vPotential,
                             XPhysicalLattice* vLattice,
                             G4int vParticleCharge){
    if(vParticleCharge < 0.){
        return fIntDensElectronNegative->GetIntegratedDensity(vPotential,
                                                              vLattice,
                                                              vParticleCharge);
    }
    else if(vParticleCharge > 0.){
        return fIntDensElectronPositive->GetIntegratedDensity(vPotential,
                                                              vLattice,
                                                              vParticleCharge);
    }
    else{
        G4cout << "XCrystalIntegratedDensityHub:: ERROR  - Charge == 0";
        G4cout << G4endl;
        return -1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::PrintOnFiles(const G4String& vFileName){
    std::string filename;
    fIntDensElectronNegative->PrintOnFile(filename=vFileName + "_neg_eld.txt");
    fIntDensElectronPositive->PrintOnFile(filename=vFileName + "_pos_eld.txt");

    fIntDensNucleiNegative->PrintOnFile(filename=vFileName + "_neg_nud.txt");
    fIntDensNucleiPositive->PrintOnFile(filename=vFileName + "_pos_nud.txt");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XCrystalIntegratedDensityHub::ReadFromFiles(const G4String& vFileName){
    std::string filename;
    fIntDensElectronNegative->ReadFromFile(filename=vFileName + "_neg_eld.txt");
    fIntDensElectronPositive->ReadFromFile(filename=vFileName + "_pos_eld.txt");
    
    fIntDensNucleiNegative->ReadFromFile(filename=vFileName + "_neg_nud.txt");
    fIntDensNucleiPositive->ReadFromFile(filename=vFileName + "_pos_nud.txt");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
