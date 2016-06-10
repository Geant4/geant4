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

#include "XVCrystalIntegratedDensity.hh"

XVCrystalIntegratedDensity::XVCrystalIntegratedDensity(){
    fNumberOfPoints = 512;
    fIntegrationPoints[0] = 32;
    fIntegrationPoints[1] = 32;
    fIntegrationPoints[2] = 32;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalIntegratedDensity::~XVCrystalIntegratedDensity(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::
SetIntegrationPoints(unsigned int vIndex,
                     unsigned int vIntegrationPoints){
    if(vIndex<3) {
        if(vIntegrationPoints > 0){
            fIntegrationPoints[vIndex] = vIntegrationPoints;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XVCrystalIntegratedDensity::
GetIntegrationPoints(unsigned int vIndex){
    if(vIndex<3) {
        return fIntegrationPoints[vIndex];
    }
    else{
        return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XVCrystalIntegratedDensity::GetIntegrationPoints(){
    return fIntegrationPoints[0]*fIntegrationPoints[1]*fIntegrationPoints[2];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::
SetNumberOfPoints(unsigned int vNumberOfPoints){
    fNumberOfPoints = vNumberOfPoints;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

unsigned int XVCrystalIntegratedDensity::GetNumberOfPoints(){
    return fNumberOfPoints;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::
SetPotential(XVCrystalCharacteristic* vPotential){
    fPotential = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XVCrystalIntegratedDensity::GetPotential(){
    return fPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::
SetDensity(XVCrystalCharacteristic* vDensity){
    fDensity = vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* XVCrystalIntegratedDensity::GetDensity(){
    return fDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::
SetXPhysicalLattice(XPhysicalLattice* vLattice){
    fLattice = vLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice*  XVCrystalIntegratedDensity::GetXPhysicalLattice(){
    return fLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::SetParticleCharge(G4int vParticleCharge){
    fParticleCharge = vParticleCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int  XVCrystalIntegratedDensity::GetParticleCharge(){
    return fParticleCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::GetStep(){
    return fPotentialRange / fNumberOfPoints;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XVCrystalIntegratedDensity::
HasBeenInitialized(XPhysicalLattice* vLattice,
                   G4int vParticleCharge){
    //now it checks only of the table is initialized,
    //it does not check if the particular crystal is initialized.
    //To be changed in the future!
    if(fTableVector->GetVectorLength() == 0) return false;
    else if(vLattice!=fLattice) return false;
    else if(vParticleCharge!=fParticleCharge) return false;
    else return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::ComputePotentialParameters(){
    fPotentialMinimum = fPotential->GetMinimum(fLattice);
    if(fParticleCharge < 0.){
        fPotentialMinimum = - fPotential->GetMaximum(fLattice);
    }
    
    fPotentialMaximum = fPotential->GetMaximum(fLattice);
    if(fParticleCharge < 0.){
        fPotentialMaximum = - fPotential->GetMinimum(fLattice);
    }
    
    fPotentialRange = std::fabs(fPotentialMaximum - fPotentialMinimum);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::InitializeTable(){
    
    ComputePotentialParameters();
    
    G4cout << "XVCrystalIntegratedDensity::InitializeTable()::";
    G4cout << "Potential Range =  " << fPotentialRange/CLHEP::eV;
    G4cout << " - Minimum = " << fPotentialMinimum / CLHEP::eV;
    G4cout << " - Maximum " << fPotentialMaximum / CLHEP::eV << G4endl;
    
    G4double vPotentialInitial = 0.;
    
    fTableVector =
    new G4PhysicsLinearVector(0.,
                               std::fabs(fPotentialMaximum - fPotentialMinimum),
                               GetNumberOfPoints());
    
    G4double vValue = 0;
    for(unsigned int i=0;i<GetNumberOfPoints();i++){
        vPotentialInitial = (fPotentialMinimum +
                             fPotentialRange * G4double(i+1)
                             / G4double(GetNumberOfPoints()));
        vValue = ComputeIntegratedDensity(vPotentialInitial,fParticleCharge);
        fTableVector->PutValue(i,vValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::
GetIntegratedDensity(G4double vPotential,
                     XPhysicalLattice* vLattice,
                     G4int vParticleCharge){
        
    G4double vPotentialModified = vPotential /
        std::fabs(G4double(vParticleCharge)) - fPotentialMinimum;
    
    // if the crystal has not been initialized return -1 -> Error!
    if(!HasBeenInitialized(vLattice,vParticleCharge))
        return -1.;
    // if the potential is higher than the maximum the average density
    // is equal to the one of the amorphous material
    else if(vPotentialModified >= std::fabs(fPotentialMaximum - fPotentialMinimum))
        return 1.;
    // if the value is less than zero (because of possible variation
    // due to centrifugal force) take the zero value
    else if(vPotentialModified < 0.) return fTableVector->Value(0.);
    //else if(vPotentialModified < 1. && vParticleCharge < 0.) return 1.;
    else{
        return fTableVector->Value(vPotentialModified);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::
ComputeIntegratedDensity(G4double vPotentialInitial,
                         G4int){
    
    unsigned int i1,i2,i3;
    i1 = i2 = i3 = 0;
    
    G4ThreeVector vPositionTemp = G4ThreeVector(0.,0.,0.);
    G4double vDensity = 0.;
    
    G4ThreeVector vSize = fLattice->GetXUnitCell()->GetSize();
    while(i1<fIntegrationPoints[2]){
        vPositionTemp.setY(G4double(G4double(i3)/
                            G4double(fIntegrationPoints[2])*vSize.z()));
        while(i1<fIntegrationPoints[1]){
            vPositionTemp.setZ(G4double(G4double(i2)/
                                G4double(fIntegrationPoints[1])*vSize.y()));
            while(i1<fIntegrationPoints[0]){
                vPositionTemp.setX(G4double(G4double(i1)/
                                G4double(fIntegrationPoints[0])*vSize.x()));
                if(fPotential->GetEC(vPositionTemp,fLattice).x()
                   < vPotentialInitial){
                    vDensity += fDensity->GetEC(vPositionTemp,fLattice).x();
                }
                i1++;
            };
            i2++;
        };
        i3++;
    };
    
    vDensity *= fLattice->GetXUnitCell()->ComputeVolume();
    vDensity /= GetIntegrationPoints();
    
    return vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalIntegratedDensity::
FindCatmullRomInterpolate(G4double &p0,
                          G4double &p1,
                          G4double &p2,
                          G4double &p3,
                          G4double &x)
{
    double a0, a1 , a2 , a3 , x2;
    
    x2 = x * x;
    a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    a1 = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
    a2 = -0.5 * p0 + 0.5 * p2;
    a3 = p1;
    
    return (a0 * x * x2 + a1 * x2 + a2 * x + a3);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::PrintOnFile(const G4String& filename){
    std::ofstream vFileOut;
    vFileOut.open(filename);
    G4double vStep = GetStep();
    
    vFileOut << "energy,dens" << std::endl;
    for(unsigned int i = 0;i < fNumberOfPoints;i++){
        vFileOut << i * vStep  / CLHEP::eV;
        vFileOut << " , ";
        vFileOut << GetIntegratedDensity((i+1) * vStep
                    + fPotentialMinimum,fLattice,fParticleCharge);
        vFileOut << std::endl;
    }
    vFileOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalIntegratedDensity::ReadFromFile(const G4String& filename){
    std::ifstream vFileIn;
    vFileIn.open(filename);
    if(!vFileIn){
        G4cout << "XVCrystalIntegratedDensity::";
        G4cout << "ReadFromFile - ERROR READING FILE!!!!! ";
        G4cout << filename << G4endl;
    }
    
    fPotentialMinimum = +DBL_MAX;
    fPotentialMaximum = -DBL_MAX;
    G4double vDensity = 0.;
    
    vFileIn >> fPotentialMinimum;
    vFileIn >> fPotentialMaximum;
    
    fPotentialMinimum *= CLHEP::eV;
    fPotentialMaximum *= CLHEP::eV;

    std::vector<G4double> fTable;
    
    while(!vFileIn.eof()){
        vFileIn >> vDensity;
        if(vDensity < 1.E-2) {
            vDensity = 1.E-2;
        }
        fTable.push_back(vDensity);
    };
    
    fNumberOfPoints = fTable.size();
    
    fTableVector =
        new G4PhysicsLinearVector(0.,
                                  std::fabs(fPotentialMaximum - fPotentialMinimum),
                                  fNumberOfPoints);
    
    for(unsigned int i=0;i<fTable.size();i++){
        fTableVector->PutValue(i,fTable.at(i));
    }
    
    fPotentialRange = std::fabs(fPotentialMaximum - fPotentialMinimum);
    
    G4cout << "XVCrystalIntegratedDensity::InitializeTable()::";
    G4cout << "Potential Range =  " << fPotentialRange/CLHEP::eV;
    G4cout << " - Minimum = " << fPotentialMinimum / CLHEP::eV;
    G4cout << " - Maximum " << fPotentialMaximum / CLHEP::eV << G4endl;
    vFileIn.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
