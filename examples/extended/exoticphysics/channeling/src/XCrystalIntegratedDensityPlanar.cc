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

#include "XCrystalIntegratedDensityPlanar.hh"

XCrystalIntegratedDensityPlanar::XCrystalIntegratedDensityPlanar(){
    SetNumberOfPoints(256);
    
    SetIntegrationPoints(0,1024);
    SetIntegrationPoints(1,1);
    SetIntegrationPoints(2,1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensityPlanar::~XCrystalIntegratedDensityPlanar(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XCrystalIntegratedDensityPlanar::ComputeIntegratedDensity(
                                                G4double vPotentialInitial,
                                                G4int vParticleCharge){
    
    unsigned int i1 = 0;
    
    G4ThreeVector vPositionTemp = G4ThreeVector(0.,0.,0.);
    G4double vInterplanarPeriod =
        GetXPhysicalLattice()->ComputeInterplanarPeriod();
    G4double vPotential = 0.;
    G4double vDensity = 0.;
    G4double xPos = 0.;
    G4double xMin = +vInterplanarPeriod;
    G4double xMax = -vInterplanarPeriod;
    
    while(i1<GetIntegrationPoints(0)){
        
        xPos = G4double(G4double(i1)/G4double(GetIntegrationPoints(0)));
        xPos *= vInterplanarPeriod;
        vPositionTemp.setX(xPos + vInterplanarPeriod / 2.);
        
        vPotential = G4double(vParticleCharge);
        vPotential *= GetPotential()->GetEC(
                        vPositionTemp,GetXPhysicalLattice()).x();
        
        if(vPotential <= vPotentialInitial){
            vDensity += GetDensity()->GetEC(vPositionTemp,
                                            GetXPhysicalLattice()).x();
            if(xPos < xMin){
                xMin = xPos;
            }
            if(xPos > xMax){
                xMax = xPos;
            }
        }
        
        i1++;
    };
    
    
    vDensity *= vInterplanarPeriod;
    vDensity *= vInterplanarPeriod/std::fabs(xMax - xMin);

    vDensity /= GetIntegrationPoints(0);
    
    return vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
