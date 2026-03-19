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

#ifndef G4CHANNELINGMATERIALDATA_HH
#define G4CHANNELINGMATERIALDATA_HH

#include "G4ios.hh"
#include "globals.hh"

#include "G4VMaterialExtension.hh"
#include "G4ChannelingECHARM.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include <unordered_map>
#include "G4PhysicsLinearVector.hh"

class G4ChannelingMaterialData : public G4VMaterialExtension
{
  public:
    
    G4ChannelingMaterialData(const G4String&);
    virtual ~G4ChannelingMaterialData();
    
    void Print() const { G4cout << "Channeling Material Data" << G4endl; }
    void SetFilename(const G4String&);
    void SetFilenameElement(const G4String&, const G4String&);
    
    G4ChannelingECHARM* GetPot() { return fPotential; }
    G4ChannelingECHARM* GetEFX() { return fElectricFieldX; }
    G4ChannelingECHARM* GetEFY() { return fElectricFieldY; }
    G4ChannelingECHARM* GetNuD() { return fNucleiDensity; }
    G4ChannelingECHARM* GetElD() { return fElectronDensity; }
        
    G4ChannelingECHARM* GetPotEl(const G4String& name) { return fPotentialElement[name]; }
    G4ChannelingECHARM* GetEFXEl(const G4String& name) { return fElectricFieldXElement[name]; }
    G4ChannelingECHARM* GetEFYEl(const G4String& name) { return fElectricFieldYElement[name]; }
    G4ChannelingECHARM* GetNuDEl(const G4String& name) { return fNucleiDensityElement[name]; }
    G4ChannelingECHARM* GetElDEl(const G4String& name) { return fElectronDensityElement[name]; }
    
    virtual G4bool IsBent() { return bIsBent; }
    
    virtual G4ThreeVector GetBR(const G4ThreeVector& v3) { return G4ThreeVector(fVectorR->Value(v3.z()),0.,0.); }
    virtual void SetBR(const G4String&);
    virtual void SetBR(G4double);

  protected:

    G4PhysicsVector* fVectorR = nullptr;
    G4bool bIsBent = false;

  private:

    G4ChannelingECHARM* fPotential = nullptr;
    G4ChannelingECHARM* fElectricFieldX = nullptr;
    G4ChannelingECHARM* fElectricFieldY = nullptr;
    G4ChannelingECHARM* fNucleiDensity = nullptr;
    G4ChannelingECHARM* fElectronDensity = nullptr;

    std::unordered_map<std::string, G4ChannelingECHARM*> fPotentialElement;
    std::unordered_map<std::string, G4ChannelingECHARM*> fElectricFieldXElement;
    std::unordered_map<std::string, G4ChannelingECHARM*> fElectricFieldYElement;
    std::unordered_map<std::string, G4ChannelingECHARM*> fNucleiDensityElement;
    std::unordered_map<std::string, G4ChannelingECHARM*> fElectronDensityElement;
};

#endif
