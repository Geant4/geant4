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
// --------------------------------------------------------------
//

#ifndef XCrystalIntegratedDensityHub_h
#define XCrystalIntegratedDensityHub_h

#include "XVCrystalCharacteristic.hh"
#include "XVCrystalIntegratedDensity.hh"
#include "G4ParticleDefinition.hh"

class XCrystalIntegratedDensityHub {

public:
    void SetDensityElectron(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetDensityElectron();

    void SetDensityNuclei(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetDensityNuclei();

    void SetPotential(XVCrystalCharacteristic*);
    XVCrystalCharacteristic* GetPotential();

    void SetXPhysicalLattice(XPhysicalLattice*);
    XPhysicalLattice* GetXPhysicalLattice();

    void SetIntegratedDensityNuclei(XVCrystalIntegratedDensity*,G4int);
    XVCrystalIntegratedDensity* GetIntegratedDensityNuclei(G4int);

    void SetIntegratedDensityElectron(XVCrystalIntegratedDensity*,G4int);
    XVCrystalIntegratedDensity* GetIntegratedDensityElectron(G4int);

    void PrintOnFiles(const G4String&);
    void ReadFromFiles(const G4String&);
    G4bool HasBeenInitialized(XPhysicalLattice*);

    G4double GetIntegratedDensityElectron(G4double,XPhysicalLattice*,G4int);
    G4double GetIntegratedDensityNuclei(G4double,XPhysicalLattice*,G4int);

public:
    void InitializeTables();
    
private:
    XPhysicalLattice* fLattice;
    XVCrystalCharacteristic* fDensityElectron;
    XVCrystalCharacteristic* fDensityNuclei;
    XVCrystalCharacteristic* fPotential;

    XVCrystalIntegratedDensity* fIntDensElectronPositive;
    XVCrystalIntegratedDensity* fIntDensNucleiPositive;

    XVCrystalIntegratedDensity* fIntDensElectronNegative;
    XVCrystalIntegratedDensity* fIntDensNucleiNegative;

public:   //Contructors
    XCrystalIntegratedDensityHub();
    ~XCrystalIntegratedDensityHub();
};

#endif
