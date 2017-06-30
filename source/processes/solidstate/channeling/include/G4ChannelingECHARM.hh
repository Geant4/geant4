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

#ifndef G4ChannelingECHARM_h
#define G4ChannelingECHARM_h

#include "G4PhysicsVector.hh"
#include "G4Physics2DVector.hh"
#include "G4ThreeVector.hh"

class G4ChannelingECHARM {

public:
    G4ChannelingECHARM(const G4String&,G4double);
    ~G4ChannelingECHARM();

    //virtual function of XVCrystalCharacteristic
    G4double GetEC(G4ThreeVector&);
    
    virtual void ReadFromECHARM(const G4String&,G4double);
    
    G4double GetMax() {return fMaximum;};
    G4double GetMin() {return fMinimum;};
    G4double GetMaxMin() {return fabs(fMaximum-fMinimum);};
    
    G4double GetIntSp(G4int index) {return fDistances[index];};
    
private:
    G4PhysicsVector* fVectorEC;
    G4Physics2DVector* fVectorEC2D;

    G4double fDistances[3];
    G4int fPoints[3];

    G4double fMaximum;
    G4double fMinimum;
};

#endif
