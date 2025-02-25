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
// Author:      Alexei Sytov

#ifndef G4CoherentPairProductionPhysics_h
#define G4CoherentPairProductionPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4CoherentPairProduction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4CoherentPairProductionPhysics : public G4VPhysicsConstructor
{
public:
    G4CoherentPairProductionPhysics(const G4String& name =
                                    "Coherent Pair Production Physics");

    ~G4CoherentPairProductionPhysics() = default;

    void ConstructParticle() override;
    void ConstructProcess() override;

    ///activate incoherent scattering
    ///(standard gamma conversion should be switched off in physics list)
    void ActivateIncoherentScattering(){fIncoherentScattering = true;}

    ///set functions

    ///set name of G4ChannelingFastSimModel from which an auto input should be performed
    void SetNameChannelingModel(const G4String& nameChannelingModel)
    {fNameChannelingModel=nameChannelingModel;}

    ///set name of G4Region to where the G4ChannelingFastSimModel is active
    void SetNameG4Region(const G4String& nameG4Region)
    {fNameRegion=nameG4Region;}

    void SetLowEnergyLimit(G4double energy){fLowEnergyLimit=energy;}
    void SetHighAngleLimit(G4double angle) {fHighAngleLimit=angle;}
    void SetPPKineticEnergyCut(G4double kineticEnergyCut) {fPPKineticEnergyCut=kineticEnergyCut;}

    /// set the number of pairs in sampling of Baier-Katkov Integral
    /// (MC integration by e+- energy and angles <=> e+- momentum)
    void SetSamplingPairsNumber(G4int nPairs){fNMCPairs = nPairs;}

    /// set the number of particle angles 1/gamma in pair production
    /// defining the width of the angular distribution of pair sampling
    /// in the Baier-Katkov Integral
    void SetChargeParticleAngleFactor(G4double chargeParticleAngleFactor)
    {fChargeParticleAngleFactor = chargeParticleAngleFactor;}

    /// set number of trajectory steps of a single particle (e- or e+)
    void SetNTrajectorySteps(G4int nTrajectorySteps)
    {fNTrajectorySteps = nTrajectorySteps;}

private:

    ///flag of simulation of incoherent scattering
    G4bool fIncoherentScattering = false;

    ///name of G4ChannelingFastSimModel from which an auto input should be performed
    G4String fNameChannelingModel = "ChannelingModel";

    ///name of G4Region to where the G4ChannelingFastSimModel is active
    G4String fNameRegion = "Crystal";

    G4double fLowEnergyLimit = 1*CLHEP::GeV;
    G4double fHighAngleLimit = 50*CLHEP::mrad;

    ///minimal kinetic energy of a charged particle produced
    G4double fPPKineticEnergyCut = 1*CLHEP::MeV;

    ///Monte Carlo statistics of e+- pair sampling in Baier-Katkov for 1 photon
    G4int fNMCPairs = 150;

    G4double fChargeParticleAngleFactor = 4; // number of particle angles 1/gamma:
        // more fChargeParticleAngleFactor => higher paramParticleAngle

    ///number of trajectory steps of a single particle (e- or e+)
    G4int fNTrajectorySteps=250;

};

#endif
