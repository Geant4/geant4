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
// Co-author:   Gianfranco Paterno (testing)
// Using the key points of G4BaierKatkov and developments of V.V. Tikhomirov,
// partially described in L. Bandiera et al. Eur. Phys. J. C 82, 699 (2022)

#ifndef G4CoherentPairProduction_h
#define G4CoherentPairProduction_h 1

#include "G4VDiscreteProcess.hh"

#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Vector/TwoVector.h>

#include "G4ChannelingFastSimCrystalData.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4CoherentPairProduction : public G4VDiscreteProcess
{
public:
    G4CoherentPairProduction(const G4String& processName = "cpp",
                             G4ProcessType aType = fElectromagnetic);

    ~G4CoherentPairProduction() = default;

    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

    G4bool IsApplicable(const G4ParticleDefinition& aPD) override
    {
        return(aPD.GetParticleName() == "gamma");
    }

    // print documentation in html format
    void ProcessDescription(std::ostream&) const override;

    ///special functions
    void Input(const G4Material* crystal,
               const G4String &lattice)
    {Input(crystal,lattice,"");}

    void Input(const G4Material* crystal,
               const G4String &lattice,
               const G4String &filePath);

    // an option to use crystal data already created outside this class
    void Input(const G4ChannelingFastSimCrystalData* crystalData);

    ///activate incoherent scattering
    ///(standard gamma conversion should be switched off in physics list)
    void ActivateIncoherentScattering(){fIncoherentScattering = true;}

    G4ChannelingFastSimCrystalData* GetCrystalData() {return fCrystalData;}

    ///get cuts
    // minimal energy for non-zero cross section
    G4double ModelMinPrimaryEnergy() { return fLowEnergyLimit;}
    G4double GetHighAngleLimit() {return fHighAngleLimit;}
    G4double GetPPKineticEnergyCut() {return fPPKineticEnergyCut;}

    /// get the number of pairs in sampling of Baier-Katkov Integral
    /// (MC integration by e+- energy and angles <=> e+- momentum)
    G4int GetSamplingPairsNumber(){return fNMCPairs;}

    /// get the number of particle angles 1/gamma in pair production
    /// defining the width of the angular distribution of pair sampling
    /// in the Baier-Katkov Integral
    G4double GetChargeParticleAngleFactor(){return fChargeParticleAngleFactor;}

    /// get number of trajectory steps of a single particle (e- or e+)
    G4double GetNTrajectorySteps(){return fNTrajectorySteps;}

    /// get effective radiation length
    /// (due to coherent process of pair production)
    /// simulated for the current photon
    G4double GetEffectiveLrad(){return fEffectiveLrad;}

    ///get the name of G4Region in which the model is applicable
    G4String GetG4RegionName() {return fG4RegionName;}

    ///set cuts
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

    ///set the name of G4Region in which the model is applicable
    void SetG4RegionName(const G4String& nameG4Region){fG4RegionName=nameG4Region;}

    G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double,
                             G4ForceCondition* condition) override;

private:

    G4int FindVectorIndex(std::vector<G4double> &myvector, G4double value);

    G4ChannelingFastSimCrystalData* fCrystalData{nullptr};

    //collection of etotal
    std::vector <CLHEP::Hep2Vector> fullVectorEtotal;

    //collection of x
    std::vector <CLHEP::Hep2Vector> fullVectorX;

    //collection of y
    std::vector <CLHEP::Hep2Vector> fullVectorY;

    //collection of tx
    std::vector <CLHEP::Hep2Vector> fullVectorTX;

    //collection of tx
    std::vector <CLHEP::Hep2Vector> fullVectorTY;

    //the vector of the discrete CDF of the production of sampling e+e- pairs
    //(in reality per distance along the photon direction)
    std::vector <G4double> fPairProductionCDFdz;

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

    ///effective radiation length (due to coherent process of pair production)
    G4double fEffectiveLrad = 0.;

    ///the name of G4Region in which the model is applicable
    G4String fG4RegionName = "Crystal";

    ///charged particle mass
    const G4double fMass = CLHEP::electron_mass_c2;

    ///flag of simulation of incoherent scattering
    G4bool fIncoherentScattering = false;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

