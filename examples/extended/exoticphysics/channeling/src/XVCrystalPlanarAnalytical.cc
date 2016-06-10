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

#include "XVCrystalPlanarAnalytical.hh"
#include "G4PhysicsLinearVector.hh"

XVCrystalPlanarAnalytical::XVCrystalPlanarAnalytical(){
    fNumberOfPlanes = 4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalPlanarAnalytical::~XVCrystalPlanarAnalytical(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::SetNumberOfPlanes(G4int vNumberOfPlanes){
    fNumberOfPlanes = vNumberOfPlanes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int XVCrystalPlanarAnalytical::GetNumberOfPlanes(){
    return fNumberOfPlanes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalPlanarAnalytical::
ComputeEC(G4ThreeVector vPositionVector,
          XPhysicalLattice* vLattice){
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    
    G4double vInterplanarDistance =
        GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod();
    
    G4double vPosition = ComputePositionInUnitCell(vPositionVector,
                                                   vLattice).x();
    
    G4double vValue = 0.;
    for(int i=-int(GetNumberOfPlanes()/2);i<=+int(GetNumberOfPlanes()/2);i++){
        vValue += ComputeECForSinglePlane( ( vPosition + G4double(i) ) *
                                          vInterplanarDistance ,
                                          vLattice );
    }
    
    return G4ThreeVector(vValue,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalPlanarAnalytical::
ComputeECFromVector(G4ThreeVector vPosition){
    G4double vInterplanarPeriod = fPhysicalLattice->ComputeInterplanarPeriod();
    G4double vX = vPosition.x();
        if (vX < 0.0) {
            vX += ((int( - vX / vInterplanarPeriod) + 1.0 )
                  * vInterplanarPeriod);
        }
        else if( vX > vInterplanarPeriod ){
            vX -= ( int( vX / vInterplanarPeriod) * vInterplanarPeriod );
        }
        return G4ThreeVector(fVectorEC->Value(vX),0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XVCrystalPlanarAnalytical::
ComputePositionInUnitCell(G4ThreeVector vPosition,
                          XPhysicalLattice* vLattice){
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();

    G4double vInterplanarPeriod =
        GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod();
    
    G4double vPositionX = vPosition.x();
    
    if((vPositionX >= 0.) &&
       (vPositionX < vInterplanarPeriod)){
        return G4ThreeVector(vPositionX/vInterplanarPeriod,0.,0.);
    }
    else if(vPositionX == vInterplanarPeriod){
        return G4ThreeVector(0.,0.,0.);
    }
    else{
        vPositionX -= std::fmod(vPosition.x(),vInterplanarPeriod)
            * vInterplanarPeriod;
        return G4ThreeVector(vPositionX/vInterplanarPeriod,0.,0.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalPlanarAnalytical::ComputeMaximum(XPhysicalLattice* vLattice){
    unsigned int vPrecision = 1024;
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vStep = GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod()
        / vPrecision;
    
    G4double vMaximum = -DBL_MAX;
    G4double vValue;
    
    for(unsigned int i=0;i<vPrecision;i++){
        if( (vValue = GetEC(G4ThreeVector(vStep * i,0.,0.),vLattice).x() )
           > vMaximum) {vMaximum = vValue;}
    }
    return vMaximum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XVCrystalPlanarAnalytical::ComputeMinimum(XPhysicalLattice* vLattice){
    unsigned int vPrecision = 1024;
    G4VPhysicalVolume* vVolume = vLattice->GetVolume();
    G4double vStep = GetXPhysicalLattice(vVolume)->ComputeInterplanarPeriod()
        / vPrecision;
    
    G4double vMinimum = +DBL_MAX;
    G4double vValue;
    
    for(unsigned int i=0;i<vPrecision;i++){
        if( (vValue = GetEC(G4ThreeVector(vStep * i,0.,0.),vLattice).x() )
           < vMinimum) {vMinimum = vValue;}
    }
    return vMinimum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::PrintOnFile(const G4String& filename,
                                            XPhysicalLattice* vLattice,
                                            G4double vUnit){
    std::ofstream vFileOut;
    vFileOut.open(filename);
    vFileOut << "pos val" << std::endl;
    
    G4int imax = 8192;
    G4double vXposition = 0.;
    G4double vInterplanarPeriod = vLattice->ComputeInterplanarPeriod();
    
    for(G4int i = 0;i<imax;i++){
        vXposition = double(i) / double(imax) * vInterplanarPeriod;
        vFileOut << vXposition / CLHEP::angstrom << " ";
        vFileOut << GetEC(G4ThreeVector(vXposition,0.,0.),vLattice).x() / vUnit;
        vFileOut << std::endl;
    }
    
    vFileOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::ReadFromFile(const G4String& filename,
                                             XPhysicalLattice*,G4double){
    std::ifstream vFileIn;
    vFileIn.open(filename);
    std::string vTempString;
    vFileIn >> vTempString >> vTempString;
 
    double vTempX;
    double vTempVal;
    std::vector<double> vTempXvector;
    std::vector<double> vTempValvector;
    
    while(!vFileIn.eof()){ 
        vFileIn >> vTempX;
        vFileIn >> vTempVal;
        vTempXvector.push_back(vTempX);
        vTempValvector.push_back(vTempVal);
    }
    
    vFileIn.close();

    G4int imax = vTempXvector.size();
    G4double vInterplanarPeriod = vTempXvector.at(imax-1);

    fVectorEC = new G4PhysicsLinearVector(0,
                                          vInterplanarPeriod*(imax-1)/imax,
                                          imax);
    for(G4int i = 0;i<imax;i++){
        fVectorEC->PutValue(i,vTempValvector.at(i));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::ReadFromECHARM(const G4String& filename,
                                               G4double vConversion){
    std::ifstream vFileIn;
    vFileIn.open(filename);
    
    G4int imax;
    G4double xmax;
    
    vFileIn >> imax;
    vFileIn >> xmax;

    xmax *= CLHEP::meter;
    fMaximum = -DBL_MAX;
    fMinimum = +DBL_MAX;

    fVectorEC = new G4PhysicsLinearVector(0,xmax,imax);
    
    for(G4int i=0;i<imax; i++){ 
        double vTempX;
        vFileIn >> vTempX;

        vTempX *= vConversion;
        if(vTempX > fMaximum) {fMaximum = vTempX;}
        if(vTempX < fMinimum) {fMinimum = vTempX;}
        fVectorEC->PutValue(i,vTempX);
    }

    G4cout << "XVCrystalPlanarAnalytical::ReadFromECHARM() - " <<
     vConversion << " " << imax << " " << xmax << " " << 
     fMinimum << " " << fMaximum << G4endl;
   
    vFileIn.close();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XVCrystalPlanarAnalytical::InitializeVector(){
    G4int imax = 4096;
    G4double vXposition = 0.;
    G4double vInterplanarPeriod = fPhysicalLattice->ComputeInterplanarPeriod();
    
    fVectorEC = new G4PhysicsLinearVector(0,
                                          vInterplanarPeriod*(imax-1)/imax,
                                          imax);
    for(G4int i = 0;i<imax;i++){
        vXposition = double(i) / double(imax) * vInterplanarPeriod;
        fVectorEC->PutValue(i,ComputeEC(G4ThreeVector(vXposition,0.,0.),
                                        fPhysicalLattice).x());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
