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

#ifndef G4Channeling_h
#define G4Channeling_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ChannelingMaterialData.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"

class G4ChannelingTrackData;

class G4Channeling : public G4VDiscreteProcess
{
public:
    
    G4Channeling();
    virtual ~G4Channeling();
    virtual G4VParticleChange* PostStepDoIt(const G4Track&,
                                            const G4Step&);
    virtual G4bool IsApplicable(const G4ParticleDefinition& aPD){
        return(aPD.GetPDGCharge() != 0.);
    };
    virtual void BuildPhysicsTable(const G4ParticleDefinition&){;};
    
protected:
    virtual G4double GetMeanFreePath(const G4Track&,
                                     G4double,
                                     G4ForceCondition*);

private:
    G4ParticleDefinition* GetParticleDefinition(const G4Track& aTrack){
        return const_cast<G4ParticleDefinition*>(aTrack.GetParticleDefinition());
    }
private:
    G4StepPoint* GetPre(const G4Track& aTrack){return aTrack.GetStep()->GetPreStepPoint();}
    G4StepPoint* GetPost(const G4Track& aTrack){return aTrack.GetStep()->GetPostStepPoint();}
    
    
private:
    G4ChannelingMaterialData* GetMatData(const G4Track& aTrack){
        G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
        if(aLV->IsExtended() == true){
            G4ExtendedMaterial* aEM = (G4ExtendedMaterial*) aTrack.GetVolume()->GetLogicalVolume()->GetMaterial();
            return (G4ChannelingMaterialData*) aEM->RetrieveExtension("channeling");
        }
        else{
            return nullptr;
        }
    }
    
    //----------------------------------------
    // Functions for the calculations of
    // parameters related to channeling
    //----------------------------------------
public:
    G4double GetCriticalAngle(const G4Track& aTrack){
        return std::sqrt(2.0*GetMatData(aTrack)->GetPot()->GetMaxMin()
                         /GetPre(aTrack)->GetTotalEnergy());}
    G4double GetOscillationPeriod(const G4Track& aTrack){
        return (CLHEP::pi * GetMatData(aTrack)->GetPot()->GetIntSp(0)
                / GetCriticalAngle(aTrack));
    }
    //----------------------------------------
    // Channeling Auxiliary Track Information
    //----------------------------------------
private:
    G4int fChannelingID;
    G4ChannelingTrackData* GetTrackData(const G4Track&);

    //----------------------------------------
    // Variables for the integration
    // of the particle trajectory
    //----------------------------------------
private:
    G4bool UpdateIntegrationStep(const G4Track&,
                                 G4ThreeVector&,
                                 G4double&);
    G4bool UpdateParameters(const G4Track&);

    void GetEF(const G4Track&,G4ThreeVector&,G4ThreeVector&);

public:
    void PosToLattice(G4StepPoint* step,G4ThreeVector&);
    
public:
    G4double GetTransverseVariationMax() {return fTransverseVariationMax;};
    void SetTransverseVariationMax(G4double aDouble) {fTransverseVariationMax = aDouble;};
    
    G4double GetTimeStepMin() {return fTimeStepMin;};
    void SetTimeStepMin(G4double aDouble) {fTimeStepMin = aDouble;};
    
private:
    G4double fTimeStepMin;
    G4double fTimeStepMax;

    G4double fTransverseVariationMax;
    
    const G4ThreeVector k010;
    G4ThreeVector fSpin;
};

#endif










