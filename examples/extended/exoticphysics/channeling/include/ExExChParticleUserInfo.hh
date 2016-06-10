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

#ifndef ChannelingUserInfo_h
#define ChannelingUserInfo_h 1

#include "globals.hh"
#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"

class ExExChParticleUserInfo : public G4VUserTrackInformation
{
    
public:
    
    ExExChParticleUserInfo();
    ~ExExChParticleUserInfo();
    
    void SetCoherentEffect(G4int flag); 
    G4int HasBeenUnderCoherentEffect();
    
    void SetNucleiDensity(G4double);
    G4double GetNucleiDensity();
    
    void SetElectronDensity(G4double);
    G4double GetElectronDensity();
    
    G4double GetNucleiDensityPreviousStep();
    G4double GetElectronDensityPreviousStep();
    void StoreDensityPreviousStep();

    G4ThreeVector GetMomentumChanneled();
    void SetMomentumChanneled(G4ThreeVector);

    G4ThreeVector GetPositionChanneled();
    void SetPositionChanneled(G4ThreeVector);

    G4double GetEnergyChanneled();
    void SetEnergyChanneled(G4double);

    G4ThreeVector GetMomentumChanneledInitial();
    void SetMomentumChanneledInitial(G4ThreeVector);
    
    G4ThreeVector GetPositionChanneledInitial();
    void SetPositionChanneledInitial(G4ThreeVector);

    G4int GetNumberOfDechanneling();
    void IncreaseNumberOfDechanneling();
    
    G4int GetInTheCrystal() {return fInTheCrystal;};
    void SetInTheCrystal(G4int aInt) {fInTheCrystal = aInt;};

private:
    
    G4int bHasBeenUnderCoherentEffect;
    //Has been in channeling in the last step

    G4double fNucleiDensity;
    //Last value of density seen by channeled particle
    G4double fNucleiDensityPreviousStep;
    
    G4double fElectronDensity;
    //Last value of density seen by channeled particle
    G4double fElectronDensityPreviousStep;

    G4ThreeVector fMomentumInChanneling;
    //Last position of the particle in the channel
    G4ThreeVector fMomentumInChannelingInitial;
    //Last position of the particle in the channel

    G4ThreeVector fPositionInChanneling;
    //Last projection fof the particle momentum in the crystal reference system
    G4ThreeVector fPositionInChannelingInitial;
    //Last projection fof the particle momentum in the crystal reference system
    
    G4int fNumberOfDechanneling;
    G4int fInTheCrystal;

};



#endif
