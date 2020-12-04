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
//

#ifndef G4ChannelingMaterialData_h
#define G4ChannelingMaterialData_h 1

#include "G4ios.hh"
#include "globals.hh"

#include "G4VMaterialExtension.hh"
#include "G4ChannelingECHARM.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include <unordered_map>
#include "G4PhysicsLinearVector.hh"

class G4ChannelingMaterialData : public G4VMaterialExtension{
public:
    
    G4ChannelingMaterialData(const G4String&);
    virtual ~G4ChannelingMaterialData();
    
public:
    void Print() const {G4cout << "Channeling Material Data" << G4endl;};
    void SetFilename(const G4String&);
    void SetFilenameElement(const G4String&,std::string);
    
public:
    G4ChannelingECHARM* GetPot() {return fPotential;};
    G4ChannelingECHARM* GetEFX() {return fElectricFieldX;};
    G4ChannelingECHARM* GetEFY() {return fElectricFieldY;};
    G4ChannelingECHARM* GetNuD() {return fNucleiDensity;};
    G4ChannelingECHARM* GetElD() {return fElectronDensity;};
        
private:
    G4ChannelingECHARM* fPotential;
    G4ChannelingECHARM* fElectricFieldX;
    G4ChannelingECHARM* fElectricFieldY;
    G4ChannelingECHARM* fNucleiDensity;
    G4ChannelingECHARM* fElectronDensity;
    
public:
    G4ChannelingECHARM* GetPotEl(std::string name) {return fPotentialElement[name];};
    G4ChannelingECHARM* GetEFXEl(std::string name) {return fElectricFieldXElement[name];};
    G4ChannelingECHARM* GetEFYEl(std::string name) {return fElectricFieldYElement[name];};
    G4ChannelingECHARM* GetNuDEl(std::string name) {return fNucleiDensityElement[name];};
    G4ChannelingECHARM* GetElDEl(std::string name) {return fElectronDensityElement[name];};
    
private:
    std::unordered_map<std::string,G4ChannelingECHARM*> fPotentialElement;
    std::unordered_map<std::string,G4ChannelingECHARM*> fElectricFieldXElement;
    std::unordered_map<std::string,G4ChannelingECHARM*> fElectricFieldYElement;
    std::unordered_map<std::string,G4ChannelingECHARM*> fNucleiDensityElement;
    std::unordered_map<std::string,G4ChannelingECHARM*> fElectronDensityElement;
    
public:
    virtual G4bool IsBent() {
        return bIsBent;
    };
    
    virtual G4ThreeVector GetBR(G4ThreeVector& v3) {
        return G4ThreeVector(fVectorR->Value(v3.z()),0.,0.);
    };
    virtual void SetBR(const G4String&);
    virtual void SetBR(G4double);

protected:
    G4PhysicsVector* fVectorR;
    G4bool bIsBent;
};

#endif










