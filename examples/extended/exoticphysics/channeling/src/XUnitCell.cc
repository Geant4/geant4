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

#include "XUnitCell.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

XUnitCell::XUnitCell(){
    fSize = G4ThreeVector(0. * CLHEP::angstrom,
                          0. * CLHEP::angstrom,
                          0. * CLHEP::angstrom);
    fAngle = G4ThreeVector(0.5 * CLHEP::pi * CLHEP::radian,
                           0.5 * CLHEP::pi * CLHEP::radian,
                           0.5 * CLHEP::pi * CLHEP::radian);
    fNumberOfBases = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XUnitCell::~XUnitCell(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XUnitCell::GetSize(){
    return fSize;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XUnitCell::GetAngle(){
    return fAngle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLogicalBase* XUnitCell::GetBase(G4int i){
    if(i<fNumberOfBases){
        return fBase[i];
    }
    else{
        G4cout << "XUnitCell::GetBase(G4int) Base "
        << i << " does not exist" << std::endl;
        return NULL;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XUnitCell::SetSize(G4ThreeVector vSize){
    fSize = vSize;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XUnitCell::SetAngle(G4ThreeVector vAngle){
    fAngle = vAngle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XUnitCell::SetBase(G4int i,XLogicalBase* base){
    if(i<fNumberOfBases){
        fBase[i] = base;
    }
    else{
        G4cout << "XUnitCell::SetBase Base(G4int,G4XLogicalBase) "
        << i << " does not exist" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XUnitCell::AddBase(XLogicalBase* base){
    fNumberOfBases++;
    //the new added basis will be in the last of the [0,fNumberOfBases-1] bases
    fBase[fNumberOfBases-1] = base;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeVolume(){
    if(IsOrthogonal()){
        return ( fSize.x()*fSize.y()*fSize.z() );
    }
    else{
        return (fSize.x()*fSize.y()*fSize.z()*std::cos(fAngle.x())*std::sin(fAngle.z()));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XUnitCell::ComputeAtomVolumeDensity(){
    G4double vAtomVolumeDensity = 0.;
    for(G4int i=0;i< fNumberOfBases;i++){
        vAtomVolumeDensity += (fBase[i]->GetElement()->GetZ()
                        *fBase[i]->GetLattice()->GetLatticeNumberOfAtoms());
    }
    vAtomVolumeDensity /= ComputeVolume();
    return vAtomVolumeDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeMillerOverSizeSquared(G4int h, G4int k, G4int l){
    return std::pow(h/fSize.x(),2.) + std::pow(k/fSize.y(),2.) + std::pow(l/fSize.z(),2.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeMillerPerSizeSquared(G4int h, G4int k, G4int l){
    return std::pow(h*fSize.x(),2.) + std::pow(k*fSize.y(),2.) + std::pow(l*fSize.z(),2.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeReciprocalVectorSquared(G4int h, G4int k, G4int l){
    G4double vReciprocalVectorSquared = 0.0;

    if(IsOrthogonal()){
        vReciprocalVectorSquared = ComputeMillerOverSizeSquared(h,k,l);
    }
    else{
        G4double vTemp[6];
        G4double vVolume = ComputeVolume();
        
        vTemp[0] = fSize.y() * fSize.z() * std::sin(fAngle.y());
        vTemp[0] = std::pow(vTemp[0] * h,2.) / vVolume;
        
        vTemp[1] = fSize.x() * fSize.z() * std::sin(fAngle.x());
        vTemp[1] = std::pow(vTemp[1] * k,2.) / vVolume;
        
        vTemp[2] = fSize.x() * fSize.y() * std::sin(fAngle.z());
        vTemp[2] = std::pow(vTemp[2] * l,2.) / vVolume;
        
        vTemp[3] = fSize.x() * fSize.y() * std::pow(fSize.z(),2.) *
            (std::cos(fAngle.x()) * std::cos(fAngle.y()) - std::cos(fAngle.z()));
        vTemp[3] *= (2. * h * k);
        
        vTemp[4] = fSize.x() * std::pow(fSize.y(),2.) * fSize.z() *
            (std::cos(fAngle.z()) * std::cos(fAngle.x()) - std::cos(fAngle.y()));
        vTemp[4] *= (2. * l * h);
        
        vTemp[5] = std::pow(fSize.x(),2.) * fSize.y() * fSize.z() *
            (std::cos(fAngle.y()) * std::cos(fAngle.z()) - std::cos(fAngle.x()));
        vTemp[5] *= (2. * k * l);

        vReciprocalVectorSquared = (vTemp[0]+vTemp[1]+vTemp[2]) / vVolume;
        vReciprocalVectorSquared += (vTemp[3]+vTemp[4]+vTemp[5]);
        vReciprocalVectorSquared /= vVolume;
    }
    
    vReciprocalVectorSquared  *= (4. * CLHEP::pi * CLHEP::pi);    
    return vReciprocalVectorSquared;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeReciprocalVector(G4int h, G4int k, G4int l){
    return std::sqrt(ComputeReciprocalVectorSquared(h,k,l));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeDirectVectorSquared(G4int h, G4int k, G4int l){
    if(IsOrthogonal()){
        return ComputeMillerPerSizeSquared(h,k,l);
    }
    else{
        double vDirectVectorSquared = 0.0;
        vDirectVectorSquared = ComputeMillerPerSizeSquared(h,k,l);
        vDirectVectorSquared += 2. * h * k * fSize.y() *
            fSize.z() * std::cos(fAngle.y()) ;
        vDirectVectorSquared += 2. * l * h * fSize.x() *
            fSize.z() * std::cos(fAngle.x()) ;
        vDirectVectorSquared += 2. * l * l * fSize.x() *
            fSize.y() * std::cos(fAngle.z()) ;
        return vDirectVectorSquared;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeDirectVector(G4int h, G4int k, G4int l){
    return std::sqrt(ComputeDirectVectorSquared(h,k,l));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeDirectPeriodSquared(G4int h, G4int k, G4int l){
    return (1./ComputeMillerOverSizeSquared(h,k,l));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double XUnitCell::ComputeDirectPeriod(G4int h, G4int k, G4int l){
    return std::sqrt(ComputeDirectPeriodSquared(h,k,l));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4bool XUnitCell::IsOrthogonal(){
    if(fAngle.x() == CLHEP::pi / 2.)
        if(fAngle.y() == CLHEP::pi / 2.)
            if(fAngle.z() == CLHEP::pi / 2.)
                return true;
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool XUnitCell::IsCubic(){
    if(IsOrthogonal())
        if(fSize.x() == fSize.y())
            if(fSize.y() == fSize.z())
                return true;
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
