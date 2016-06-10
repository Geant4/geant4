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

#ifndef ExExChProcessChanneling_h
#define ExExChProcessChanneling_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

#include "XLatticeManager3.hh"
#include "XVCrystalCharacteristic.hh"
#include "XCrystalIntegratedDensityHub.hh"

#include "G4VPhysicalVolume.hh"
#include "XPhysicalLattice.hh"
#include "ExExChParticleUserInfo.hh"

class ExExChProcessChanneling : public G4VDiscreteProcess
{
public:
    
    ExExChProcessChanneling(const G4String& processName =
                                  "ExExChProcessChanneling" );
    
    virtual ~ExExChProcessChanneling();
    
    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
        
    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    
    virtual void BuildPhysicsTable(const G4ParticleDefinition&);
    
protected:
    
    virtual G4double GetMeanFreePath(const G4Track&,
                                     G4double,
                                     G4ForceCondition* );
    
private:
    
    G4double GetChannelingMeanFreePath(const G4Track&);
    
public:
    XVCrystalCharacteristic* GetPotential();
    void SetPotential(XVCrystalCharacteristic*);

    XVCrystalCharacteristic* GetElectricField();
    void SetElectricField(XVCrystalCharacteristic*);

    XCrystalIntegratedDensityHub* GetIntegratedDensity();
    void SetIntegratedDensity(XCrystalIntegratedDensityHub*);

    XVCrystalCharacteristic* GetNucleiDensity();
    void SetNucleiDensity(XVCrystalCharacteristic*);

    XVCrystalCharacteristic* GetElectronDensity();
    void SetElectronDensity(XVCrystalCharacteristic*);

    G4double GetTransverseVariationMax(){
        return fTransverseVariationMax;
        };
        
    void SetTransverseVariationMax(G4double aDouble){
        fTransverseVariationMax = aDouble;
        };

    G4double GetTimeStepMin() {return fTimeStepMin;};
    void SetTimeStepMin(G4double aDouble) {fTimeStepMin = aDouble;};

    void ReadFromFileCharacteristics(G4bool);
        
    void SetFileCharacteristicsName(const G4String& vFilename)
    {fFileCharacteristicsName = vFilename;};
    G4bool GetFileCharacteristicsName() {return fFileCharacteristicsName;};
    
private:
    void UpdateParameters(const G4Track&);
    G4bool UpdateInitialParameters(const G4Track&);
    void ResetDensity(const G4Track&);

    G4double ComputeCriticalEnergyBent(const G4Track&);
    G4double ComputeCriticalEnergyMinimumBent(const G4Track&);
    G4double ComputePotentialEnergyBent(const G4Track&);
    G4ThreeVector ComputeTransverseEnergyBent(const G4Track&);
    
    G4ThreeVector ComputePositionInTheCrystal(G4StepPoint*,const G4Track&);
    G4StepPoint* CheckStepPointLatticeForVolume(G4StepPoint*,const G4Track&);
    G4StepPoint* CheckStepPointLatticeForPosition(G4StepPoint*,const G4Track&);
    
    G4ThreeVector ComputeTransverseEnergy(const G4Track&);
    G4ThreeVector ComputeKineticEnergy(const G4Track&);
    G4ThreeVector ComputePotentialEnergy(const G4Track&);
    G4ThreeVector ComputeMomentum(const G4Track&,G4StepPoint*);
    G4ThreeVector ComputeCentrifugalEnergy(const G4Track&,G4ThreeVector);
    G4ThreeVector ComputeCentrifugalEnergyMaximumVariation(const G4Track&);

    G4double ComputeCriticalEnergyMaximum(const G4Track&);
    G4double ComputeCriticalEnergyMinimum(const G4Track&);
    G4double ComputeCriticalAngle(const G4Track&);
    G4double ComputeCriticalRadius(const G4Track&);
    G4double ComputeOscillationPeriod(const G4Track&);
    G4ThreeVector ComputePotentialWellCentre(const G4Track&);
    
    G4bool HasLattice(const G4Track&);
    G4bool HasLatticeOnBoundary(const G4Track&);
    G4bool HasLatticeOnBoundaryPre(const G4Track&);
    G4bool HasLatticeOnBoundaryPost(const G4Track&);

    G4bool ParticleIsNegative(const G4Track&);
    G4bool ParticleIsNotOnBoundaryPre(const G4Track&);
    G4bool ParticleIsNotOnBoundaryPost(const G4Track&);
    G4bool ParticleIsNotOnBoundary(const G4Track&);

    void ComputeCrystalCharacteristic(const G4Track&);

private:
    //binding methods
    XPhysicalLattice* GetXPL(const G4Track&);
    G4VPhysicalVolume* GetVolume(const G4Track&);
    ExExChParticleUserInfo* GetInfo(const G4Track&);
    G4ParticleDefinition* GetParticleDefinition(const G4Track& aTrack);
    
private:
    // hide assignment operator as private
    ExExChProcessChanneling(ExExChProcessChanneling&);
    ExExChProcessChanneling& operator=(
                            const ExExChProcessChanneling& right);
    
private:
    XLatticeManager3* fLatticeManager;
      
    XVCrystalCharacteristic* fPotentialEnergy;
    XVCrystalCharacteristic* fElectricField;
    XCrystalIntegratedDensityHub* fIntegratedDensity;

    XVCrystalCharacteristic* fNucleiDensity;
    XVCrystalCharacteristic* fElectronDensity;

    G4String fFileCharacteristicsName;
    
private:
    G4bool UpdateIntegrationStep(const G4Track&,G4ThreeVector&);
    G4double fTimeStep;
    G4double fTimeStepMin;
    G4double fTimeStepMax;
    G4double fTimeStepTotal;

    G4bool bHasToComputeTrajectory;
    G4double bPointYPost;
    G4double bPointYPre;
    
    G4double fTransverseVariationMax;
    G4double fIntegrationPeriod;
};

#endif










