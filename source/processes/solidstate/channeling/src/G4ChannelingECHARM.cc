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

#include "G4ChannelingECHARM.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4Physics2DVector.hh"
#include "G4SystemOfUnits.hh"

G4ChannelingECHARM::G4ChannelingECHARM(const G4String& fileName,G4double vConversion):
fVectorEC(0),
fDistances{0.,0.,0.},
fPoints{0,0,0},
fMaximum(-DBL_MAX),
fMinimum(DBL_MAX){
    fDistances[0] = 0;
    fDistances[1] = 0;
    fDistances[2] = 0;
    fPoints[0] = 0;
    fPoints[1] = 0;
    fPoints[2] = 0;
    ReadFromECHARM(fileName,vConversion);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ChannelingECHARM::~G4ChannelingECHARM(){
     delete(fVectorEC);
     delete(fVectorEC2D);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ChannelingECHARM::GetEC(G4ThreeVector& vPosition){
    G4double vX = vPosition.x();
    if (vX < 0.0) {
        vX += ((int( - vX / fDistances[0]) + 1.0 ) * fDistances[0]);
    }
    else if( vX > fDistances[0] ){
        vX -= ( int( vX / fDistances[0]) * fDistances[0] );
    }
    if(fPoints[1]==1){
        return fVectorEC->Value(vX);
    }
    else{
        G4double vY = vPosition.y();
        if (vY < 0.0) {
            vY += ((int( - vY / fDistances[1]) + 1.0 ) * fDistances[1]);
        }
        else if( vY > fDistances[1] ){
            vY -= ( int( vY / fDistances[1]) * fDistances[1] );
        }
        return fVectorEC2D->Value((vX),(vY));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingECHARM::ReadFromECHARM(const G4String& filename,
                                        G4double vConversion){
    std::ifstream vFileIn;
    vFileIn.open(filename);
    
    vFileIn >> fPoints[0] >> fPoints[1] >> fPoints[2];
    vFileIn >> fDistances[0] >> fDistances[1] >> fDistances[2];

    fDistances[0] *= CLHEP::meter;
    fDistances[1] *= CLHEP::meter;
    fDistances[2] *= CLHEP::meter;
    fMaximum = -DBL_MAX;
    fMinimum = +DBL_MAX;

    if(fPoints[1]<1){
        G4ExceptionDescription ed;
        ed << "No Points not found !" << G4endl;
        G4Exception("G4ChannelingECHARM::ReadFromECHARM(...)",
                  "G4ChannelingECHARM",
                  FatalException,
                  ed);
        return;
    }
    else if(fPoints[1]==1){
        fVectorEC = new G4PhysicsLinearVector(0,fDistances[0],fPoints[0]);
    }
    else{
        fVectorEC2D = new G4Physics2DVector(fPoints[0],fPoints[1]);
    }
    G4double stepX = fDistances[0]/fPoints[0];
    G4double stepY = fDistances[1]/fPoints[1];
    for(G4int i1=0;i1<fPoints[1]; i1++){
        if(fPoints[1]>1){
            fVectorEC2D->PutY(i1,i1*stepY);
        }
        for(G4int i0=0;i0<fPoints[0]; i0++){
        double vTempX;
        vFileIn >> vTempX;

        vTempX *= vConversion;
        if(vTempX > fMaximum) {fMaximum = vTempX;}
        if(vTempX < fMinimum) {fMinimum = vTempX;}
        if(fPoints[1]==1){
            fVectorEC->PutValue(i0,vTempX);
        }
        else{
            fVectorEC2D->PutValue(i0,i1,vTempX);
            fVectorEC2D->PutX(i0,i0*stepX);
        }
    }
    }
    G4cout << "G4ChannelingECHARM::ReadFromECHARM() - " << vConversion << " " << fPoints[0] << " " << fDistances[0] << " " << fPoints[1] << " " << fDistances[1] << " " << fMinimum << " " << fMaximum << G4endl;

    vFileIn.close();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
