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

#include "ExExChProcessChanneling.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4RandomTools.hh"


#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "XLatticeManager3.hh"

#include "XLogicalAtomicLattice.hh"
#include "XLogicalAtomicLatticeDiamond.hh"
#include "XLogicalBase.hh"
#include "XUnitCell.hh"

#include "ExExChParticleUserInfo.hh"

ExExChProcessChanneling::
ExExChProcessChanneling(const G4String& aName):
G4VDiscreteProcess(aName){
    fLatticeManager = XLatticeManager3::GetXLatticeManager();
    
    G4double kCarTolerance =
        G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    G4cout<<"\n ExExChProcessChanneling::Constructor:";
    G4cout<<"Geometry surface tolerance is: ";
    G4cout<< kCarTolerance / mm << " mm"<<G4endl;
    if(verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    
    fFileCharacteristicsName = "";
    fTimeStepMax = 0.;
    fTimeStepMin = 2.E2 * CLHEP::angstrom;
    fTransverseVariationMax = 2.E-2 * CLHEP::angstrom;

    bPointYPre = -1.;
    bPointYPost = -1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExExChProcessChanneling::~ExExChProcessChanneling(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExExChProcessChanneling::
ExExChProcessChanneling(ExExChProcessChanneling& right):
G4VDiscreteProcess(right){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ExExChProcessChanneling::GetPotential(){
    return fPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::
SetPotential(XVCrystalCharacteristic* vPotential){
    fPotentialEnergy = vPotential;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ExExChProcessChanneling::GetElectricField(){
    return fElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::
SetElectricField(XVCrystalCharacteristic* vElectricField){
    fElectricField = vElectricField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XCrystalIntegratedDensityHub* ExExChProcessChanneling::
GetIntegratedDensity(){
    return fIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::
SetIntegratedDensity(XCrystalIntegratedDensityHub* vIntegratedDensity){
    fIntegratedDensity = vIntegratedDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ExExChProcessChanneling::GetNucleiDensity(){
    return fNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::
SetNucleiDensity(XVCrystalCharacteristic* vNucleiDensity){
    fNucleiDensity = vNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XVCrystalCharacteristic* ExExChProcessChanneling::GetElectronDensity(){
    return fElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::
SetElectronDensity(XVCrystalCharacteristic* vElectronDensity){
    fElectronDensity = vElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::UpdateInitialParameters(const G4Track& aTrack){
    
    if(GetInfo(aTrack)->GetPositionChanneled().x() == DBL_MAX){
        // when the particle enters the crystal the position in the channel
        //is randomly generated using a uniform distribution
        G4double vXposition = G4UniformRand() *
        GetXPL(aTrack)->ComputeInterplanarPeriod();
        
        //vXposition = 1.0 * CLHEP::angstrom;
    
        //initial position in the channel is stored
        GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(vXposition,
                                                            0.,
                                                            0.));
        GetInfo(aTrack)->SetPositionChanneledInitial(G4ThreeVector(vXposition,
                                                                   0.,
                                                                   0.));
    }
    
    if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
        // the first time it enter the crystal we take the momentum
        // for the post step which is the only one in the crystal
        G4ThreeVector vMomentum =
        ComputeMomentum(aTrack,aTrack.GetStep()->GetPostStepPoint());
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
        GetInfo(aTrack)->SetMomentumChanneledInitial(vMomentum);
        return true;
    }
    
    return false;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::UpdateParameters(const G4Track& aTrack){
    
    if(fIntegratedDensity->HasBeenInitialized(GetXPL(aTrack))
       == false){
        ComputeCrystalCharacteristic(aTrack);
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::";
        G4cout<<"fIntegratedDensity->Initialized" << G4endl;
    }
    
    if(UpdateInitialParameters(aTrack) == false){
        G4ThreeVector vMomentumNew = 
                ComputeMomentum(aTrack,
                                aTrack.GetStep()->GetPreStepPoint());
        GetInfo(aTrack)->SetMomentumChanneled(vMomentumNew);
    }
    G4ThreeVector vPositionPost =
    ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),aTrack);
    G4ThreeVector vPositionPre =
    ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
    
    bHasToComputeTrajectory = true;

    if(vPositionPost.y() == bPointYPost &&
       vPositionPre.y() == bPointYPre){
        bHasToComputeTrajectory = false;
    }
    else{
        if(GetXPL(aTrack)->IsBent()){
                bPointYPre = vPositionPre.y();
                bPointYPost = vPositionPost.y();
        }
        else{
                bPointYPre = vPositionPre.z();
                bPointYPost = vPositionPost.z();
        }
        }

    G4ThreeVector fMom = GetInfo(aTrack)->GetMomentumChanneled();
    G4ThreeVector fPos = GetInfo(aTrack)->GetPositionChanneled();
    G4ThreeVector fMomHalf = GetInfo(aTrack)->GetMomentumChanneled();
    G4ThreeVector fPosHalf = GetInfo(aTrack)->GetPositionChanneled();

    if(GetXPL(aTrack)->IsBent()){
            fIntegrationPeriod = (vPositionPost.phi() - vPositionPre.phi())*
                    GetXPL(aTrack)->GetCurvatureRadius().x();
            fIntegrationPeriod = vPositionPost.y() - vPositionPre.y();
    }
    else{
            fIntegrationPeriod = vPositionPost.z() - vPositionPre.z();
    }

    fTimeStepTotal = 0.;

    G4double vNucleiDensity=0.;
    G4double vElectronDensity=0.;

    if(fIntegrationPeriod>0. && bHasToComputeTrajectory==true){
            G4double kBeta = 0.;
            G4double kPos = 0.;
            G4double kMom = 0.;
            G4double kBR = 0.;
            G4double Z = 0.;
        do{
            UpdateIntegrationStep(aTrack,fMom);

            fPosHalf = fPos;
            fMomHalf = fMom;

            kBeta = aTrack.GetVelocity()/c_light;
            if(fMom.z()!=0.){
              kPos = fTimeStep / fMom.z();
            }
            else{
              kPos = fTimeStep / 1.E-20;
            }
            kMom = fTimeStep / kBeta;
            kBR = fTimeStep * (fMom.z() * kBeta);;
            Z = GetParticleDefinition(aTrack)->GetPDGCharge();

            fPosHalf += (fMom * kPos * 0.5);
            fMomHalf += 
                    (GetElectricField()->GetEC(fPos,GetXPL(aTrack))
                    * Z  * kMom * 0.5);

            if(GetXPL(aTrack)->IsBent()){
                    G4double temp =
                            fMomHalf.x() + kBR * 0.5 /
                            (GetXPL(aTrack)->GetCurvatureRadius()).x();
                fMomHalf.setX(temp);
            }

            fPos += (fMomHalf * kPos);
            fMom += 
                    (GetElectricField()->GetEC(fPosHalf,GetXPL(aTrack))
                     * Z * kMom );

            if(GetXPL(aTrack)->IsBent()){
                    G4double temp = 
                            fMom.x() + kBR / 
                            (GetXPL(aTrack)->GetCurvatureRadius()).x();
               fMom.setX(temp);
            }

            fTimeStepTotal += fTimeStep;

            vNucleiDensity += 
                 (fTimeStep * 
                 (GetNucleiDensity()->GetEC(fPos,GetXPL(aTrack)).x()
                 +GetNucleiDensity()->GetEC(fPos,GetXPL(aTrack)).x()
                 ) * 0.5);
            vElectronDensity +=
                 (fTimeStep * (
                 GetElectronDensity()->GetEC(fPos,GetXPL(aTrack)).x() + 
                 GetElectronDensity()->GetEC(fPos,GetXPL(aTrack)).x()) * 0.5);


        } while(fTimeStepTotal<fIntegrationPeriod);

        vNucleiDensity /= fIntegrationPeriod;
        vElectronDensity /= fIntegrationPeriod;
        GetInfo(aTrack)->SetNucleiDensity(vNucleiDensity);
        GetInfo(aTrack)->SetElectronDensity(vElectronDensity);
    }
    else{
        ResetDensity(aTrack);
    }
    
    GetInfo(aTrack)->SetMomentumChanneled(fMom);
    GetInfo(aTrack)->SetPositionChanneled(fPos);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::ResetDensity(const G4Track& aTrack){
    GetInfo(aTrack)->SetNucleiDensity(1.);
    GetInfo(aTrack)->SetElectronDensity(1.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4ThreeVector ExExChProcessChanneling::
ComputePositionInTheCrystal(G4StepPoint* vStep,const G4Track& aTrack){
    
    G4StepPoint* vStepVol = CheckStepPointLatticeForVolume(vStep,aTrack);
    G4StepPoint* vStepPos = CheckStepPointLatticeForPosition(vStep,aTrack);
    
    G4TouchableHistory* theTouchable =
        (G4TouchableHistory*)(vStepVol->GetTouchable());
    G4ThreeVector vWorldPos = vStepPos->GetPosition();
    G4ThreeVector vLocalPos =
        theTouchable->GetHistory()->GetTopTransform().TransformPoint(vWorldPos);
    
    
    if(GetXPL(aTrack)->IsBent() == false){
        G4Box* vXtalSolid =
            (G4Box*) vStepVol->GetPhysicalVolume()
            ->GetLogicalVolume()->GetSolid();
        vLocalPos += G4ThreeVector(vXtalSolid->GetXHalfLength(),
                                   vXtalSolid->GetYHalfLength(),
                                   vXtalSolid->GetZHalfLength());
    }
    else{
        G4Tubs* vXtalSolid =
            (G4Tubs*) vStepVol->GetPhysicalVolume()
            ->GetLogicalVolume()->GetSolid();
        vLocalPos.rotateZ(-vXtalSolid->GetStartPhiAngle());
    }
    return vLocalPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4StepPoint* ExExChProcessChanneling::
CheckStepPointLatticeForVolume(G4StepPoint* vStep, const G4Track& aTrack){
    
    G4StepPoint* vStepPre = aTrack.GetStep()->GetPreStepPoint();
    G4StepPoint* vStepPost = aTrack.GetStep()->GetPostStepPoint();
    
    if( fLatticeManager->HasLattice(vStep->GetPhysicalVolume()) ) {
        return vStep;
    }
    else if(fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume())
            == false
            && fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume())
            == true
            && vStep == vStepPost &&
            vStepPost->GetStepStatus() == fGeomBoundary) {
        return vStepPre;
    }
    else if(fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume())
            == false &&
            fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume())
            == true &&
            vStep == vStepPre &&
            vStepPre->GetStepStatus() == fGeomBoundary) {
        return vStepPost;
    }
    else{
        return vStep;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4StepPoint* ExExChProcessChanneling::
CheckStepPointLatticeForPosition(G4StepPoint* vStep, const G4Track& aTrack){
    
    G4StepPoint* vStepPre = aTrack.GetStep()->GetPreStepPoint();
    G4StepPoint* vStepPost = aTrack.GetStep()->GetPostStepPoint();
    
    if( fLatticeManager->HasLattice(vStep->GetPhysicalVolume()) ) {
        return vStep;
    }
    else if(fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume())
            == false &&
            fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume())
            == true &&
            vStep == vStepPost &&
            vStepPost->GetStepStatus() == fGeomBoundary) {
        return vStepPost;
    }
    else if(fLatticeManager->HasLattice(vStepPre->GetPhysicalVolume())
            == false &&
            fLatticeManager->HasLattice(vStepPost->GetPhysicalVolume())
            == true &&
            vStep == vStepPre &&
            vStepPre->GetStepStatus() == fGeomBoundary) {
        return vStepPre;
    }
    else{
        return vStep;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
UpdateIntegrationStep(const G4Track& aTrack,G4ThreeVector& vMom){
    
    if(vMom.x() != 0.0 || vMom.y() != 0.0){
        double xy2 = vMom.x() * vMom.x() + vMom.y()*vMom.y();
        
        if(xy2!=0.){
            fTimeStep =
                    std::fabs(fTransverseVariationMax * 
                    aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy() /
                    std::pow(xy2,0.5));
            if(fTimeStep < fTimeStepMin) fTimeStep = fTimeStepMin;
            else{
                fTimeStepMax = std::sqrt( fTransverseVariationMax *
                aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy() /
                std::fabs(fElectricField->GetMaximum(GetXPL(aTrack))));
                
                if(fTimeStep > fTimeStepMax) fTimeStep = fTimeStepMax;
            }
        }
        else{
            fTimeStep = fTimeStepMin;
        }
        
        if(fTimeStep + fTimeStepTotal > fIntegrationPeriod){
            fTimeStep = fIntegrationPeriod - fTimeStepTotal;
        }
        
        return true;
    }
    else{
        fTimeStep = fTimeStepMin;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ExExChProcessChanneling::
GetChannelingMeanFreePath(const G4Track& aTrack){
    //----------------------------------------
    // return the channeling MFP
    //----------------------------------------
    
    
    G4double vMFP = 0.1 * ComputeOscillationPeriod(aTrack);
    
    if(HasLatticeOnBoundaryPre(aTrack) == true){
        vMFP = 0.001 * ComputeOscillationPeriod(aTrack);
    }
    
    return vMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ExExChProcessChanneling::
GetMeanFreePath(const G4Track& aTrack,
                G4double, // previousStepSize
                G4ForceCondition* condition){
    
    //----------------------------------------
    // the condition is forced to check if
    // the volume has a lattice at each step.
    // if it hasn't, return DBL_MAX
    //----------------------------------------
    
    *condition = Forced;
    
    if(HasLattice(aTrack)){
        GetInfo(aTrack)->SetInTheCrystal(true);
        return GetChannelingMeanFreePath(aTrack);
    }
    else{
        GetInfo(aTrack)->SetInTheCrystal(false);
        return DBL_MAX;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* ExExChProcessChanneling::
PostStepDoIt(const G4Track& aTrack,
             const G4Step&){
    
    //----------------------------------------
    // check if the volume has a lattice
    // and if the particle is in channeling.
    // If it is so, the particle is forced
    // to follow the channeling plane
    // direction. If the particle has
    // dechanneled or exited the crystal,
    // the outgoing angle is evaluated
    //----------------------------------------
    
    aParticleChange.Initialize(aTrack);
    
    GetInfo(aTrack)->StoreDensityPreviousStep();
    
    if((HasLattice(aTrack) == true) &&
            (HasLatticeOnBoundaryPost(aTrack) == false)){
        UpdateParameters(aTrack);
        
        G4ThreeVector vMomentum = 
                GetInfo(aTrack)->GetMomentumChanneled().unit();

        G4ThreeVector vPosition;
        vPosition =
             ComputePositionInTheCrystal(aTrack.GetStep()->GetPostStepPoint(),
             aTrack);

        GetXPL(aTrack)->ProjectMomentumVectorFromLatticeToWorld(vMomentum,
                                                                vPosition);

        aParticleChange.ProposeMomentumDirection(vMomentum.unit());
    }
    else{
        // if the volume has no lattice it resets the density factors
        ResetDensity(aTrack);
        GetInfo(aTrack)->SetMomentumChanneled(G4ThreeVector(DBL_MAX,
                                                            DBL_MAX,
                                                            DBL_MAX));
                                                            
        GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(DBL_MAX,
                                                            DBL_MAX,
                                                            DBL_MAX));
    }
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ExExChProcessChanneling::
ComputeTransverseEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vTransverseEnergy = ComputePotentialEnergy(aTrack)
        + ComputeKineticEnergy(aTrack);
    return vTransverseEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4ThreeVector ExExChProcessChanneling::
ComputeKineticEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle kinetic energy
    // in the crystal reference system
    //----------------------------------------
    
    G4double vTotalEnergy =
        aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    
    G4ThreeVector vMom = GetInfo(aTrack)->GetMomentumChanneled();
    
    G4ThreeVector vKineticEnergy = 0.5 *
        G4ThreeVector((vMom.x() * vMom.x()) / vTotalEnergy,
                      0.,
                      (vMom.z() * vMom.z())  / vTotalEnergy);
    
    return vKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ExExChProcessChanneling::
ComputePotentialEnergy(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle transverse energy
    // in the crystal reference system
    //----------------------------------------
    
    
    G4ThreeVector vPotentialEnergy =
        fPotentialEnergy->GetEC(GetInfo(aTrack)->GetPositionChanneled(),
                                GetXPL(aTrack));
    
    vPotentialEnergy *= GetParticleDefinition(aTrack)->GetPDGCharge();
    
    return vPotentialEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ExExChProcessChanneling::
ComputeCentrifugalEnergy(const G4Track& aTrack,G4ThreeVector vPosition){
    //----------------------------------------
    // compute the transverse energy variation
    // in the crystal reference system
    //----------------------------------------
    
    G4double vTotalEnergy =
        aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4double vPositionX = vPosition.x();
    
    if(ParticleIsNegative(aTrack) && false){
        vPositionX -=
        GetXPL(aTrack)->ComputeInterplanarPeriod() * 0.5;
    }
    G4ThreeVector vEnergyVariation = G4ThreeVector();;
    if(GetXPL(aTrack)->IsBent()){ 
    vEnergyVariation =
        G4ThreeVector(vTotalEnergy * vPositionX /
                      GetXPL(aTrack)->GetCurvatureRadius().x(),
                      0.,
                      0.);
    }
    
    return vEnergyVariation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ExExChProcessChanneling::
ComputeMomentum(const G4Track& aTrack,G4StepPoint* vStep){
    //----------------------------------------
    // compute the particle momentum
    // in the crystal reference system
    //----------------------------------------
    
    G4ThreeVector vPosition = ComputePositionInTheCrystal(vStep,aTrack);
    
    G4ThreeVector vMomentum = aTrack.GetMomentum();

    GetXPL(aTrack)->
    ProjectMomentumVectorFromWorldToLattice(vMomentum,vPosition);

    return vMomentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ExExChProcessChanneling::
ComputeCentrifugalEnergyMaximumVariation(const G4Track& aTrack){
    //----------------------------------------
    // compute the transverse energy variation
    // in the crystal reference system
    //----------------------------------------
    
    G4double vTotalEnergy =
        aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4ThreeVector vEnergyVariation = G4ThreeVector();
    if(GetXPL(aTrack)->IsBent()){
    vEnergyVariation = G4ThreeVector(vTotalEnergy *
                    GetXPL(aTrack)->ComputeInterplanarPeriod() /
                    GetXPL(aTrack)->GetCurvatureRadius().x(),
                    0.,
                    0.);
    }
    
    return vEnergyVariation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double ExExChProcessChanneling::
ComputeCriticalEnergyMaximum(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy
    // for channeling
    //----------------------------------------
    
    G4double vCriticalEnergy = 0.;
    
    if(ParticleIsNegative(aTrack)){
        vCriticalEnergy =
            - fPotentialEnergy->GetMinimum(GetXPL(aTrack));
    }
    else{
        vCriticalEnergy =
            + fPotentialEnergy->GetMaximum(GetXPL(aTrack));
    }
    
    vCriticalEnergy *= std::fabs(GetParticleDefinition(aTrack)->GetPDGCharge());
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ExExChProcessChanneling::
ComputeCriticalEnergyMinimum(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical energy minimum
    // for channeling
    //----------------------------------------
    
    G4double vCriticalEnergy = 0.;
    
    if(ParticleIsNegative(aTrack)){
        vCriticalEnergy =
            - fPotentialEnergy->GetMaximum(GetXPL(aTrack));
    }
    else{
        vCriticalEnergy =
            + fPotentialEnergy->GetMinimum(GetXPL(aTrack));
    }
    
    vCriticalEnergy *= std::fabs(GetParticleDefinition(aTrack)->GetPDGCharge());
    
    return vCriticalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ExExChProcessChanneling::
ComputeCriticalAngle(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical angle
    // for chenneling
    //----------------------------------------
    
    G4double vTotalEnergy =
        aTrack.GetStep()->GetPostStepPoint()->GetTotalEnergy();
    G4double vCriticalAngle = std::pow( 2.0 *
                    std::fabs( ( ComputeCriticalEnergyMaximum(aTrack)
                           - ComputeCriticalEnergyMinimum(aTrack) )
                           / vTotalEnergy ) , 0.5);
    return vCriticalAngle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ExExChProcessChanneling::
ComputeOscillationPeriod(const G4Track& aTrack){
    //----------------------------------------
    // compute the particle oscillation
    // period in the crystal channel
    //----------------------------------------
    
    G4double vInterplanarPeriod =
        GetXPL(aTrack)->ComputeInterplanarPeriod();
    G4double vOscillationPeriod =
        CLHEP::pi * vInterplanarPeriod / ComputeCriticalAngle(aTrack);
    return vOscillationPeriod;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double ExExChProcessChanneling::
ComputeCriticalRadius(const G4Track& aTrack){
    //----------------------------------------
    // compute the critical radius
    // for channeling
    //----------------------------------------
    
    G4double vTotalEnergy =
        aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
    
    G4double vCriticalRadius = 1.E-20;
    if(fElectricField->GetMaximum(GetXPL(aTrack)) != 0.){
        vCriticalRadius = vTotalEnergy 
                        / fElectricField->GetMaximum(GetXPL(aTrack));
    }
    return vCriticalRadius;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ExExChProcessChanneling::
ComputePotentialWellCentre(const G4Track& aTrack){
    //----------------------------------------
    // compute the central point of
    // the potential well for channeling
    //----------------------------------------
    
    G4double vInterplanarPeriodHalf =
        0.5 * GetXPL(aTrack)->ComputeInterplanarPeriod();
    
    G4double vCentreX = vInterplanarPeriodHalf;
    
    if(GetXPL(aTrack)->IsBent()){
        G4double vTotalEnergy =
            aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
        
        G4double vPotentialWellDepth =
            ComputeCriticalEnergyMaximum(aTrack)
            - (ComputeCriticalEnergyMinimum(aTrack));
        
        vCentreX *= (1. - 0.5 * vTotalEnergy /
                     vPotentialWellDepth /
                     GetXPL(aTrack)->GetCurvatureRadius().x() *
                     vInterplanarPeriodHalf );
    }
    
    G4ThreeVector vCentre = G4ThreeVector(vCentreX,0.,0.);
    
    
    return vCentre;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool ExExChProcessChanneling::
IsApplicable(const G4ParticleDefinition& aPD){
    return(aPD.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::
BuildPhysicsTable(const G4ParticleDefinition&){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* ExExChProcessChanneling::
GetXPL(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(
                aTrack.GetStep()->GetPostStepPoint()->GetPhysicalVolume())){
        return fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()
                            ->GetPostStepPoint()->GetPhysicalVolume());
    }
    else if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPreStepPoint()
                    ->GetPhysicalVolume()) &&
            aTrack.GetStep()->GetPostStepPoint()->GetStepStatus()
            == fGeomBoundary) {
        return fLatticeManager->GetXPhysicalLattice(aTrack.GetStep()
                    ->GetPreStepPoint()->GetPhysicalVolume());
    }
    else{
        G4cout << "LATTICE NOT FOUND: ERROR" << G4endl;
        return NULL;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::HasLattice(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()
                    ->GetPostStepPoint()->GetPhysicalVolume())){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
HasLatticeOnBoundary(const G4Track& aTrack){
    if(HasLatticeOnBoundaryPost(aTrack) || HasLatticeOnBoundaryPre(aTrack)) {
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
HasLatticeOnBoundaryPre(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()->GetPreStepPoint()->
                GetPhysicalVolume()) &&
       aTrack.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
HasLatticeOnBoundaryPost(const G4Track& aTrack){
    if(fLatticeManager->HasLattice(aTrack.GetStep()->
                    GetPostStepPoint()->GetPhysicalVolume()) &&
       aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::ParticleIsNegative(const G4Track& aTrack){
    if(GetParticleDefinition(aTrack)->GetPDGCharge() < 0.) {
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
ParticleIsNotOnBoundaryPre(const G4Track& aTrack){
    if(aTrack.GetStep()->GetPreStepPoint()->GetStepStatus() != fGeomBoundary){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
ParticleIsNotOnBoundaryPost(const G4Track& aTrack){
    if(aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() != fGeomBoundary){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
ParticleIsNotOnBoundary(const G4Track& aTrack){
    if(ParticleIsNotOnBoundaryPost(aTrack) ||
       ParticleIsNotOnBoundaryPre(aTrack)){
        return true;
    }
    else{
        return false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExExChParticleUserInfo* ExExChProcessChanneling::
GetInfo(const G4Track& aTrack){
    return (ExExChParticleUserInfo*) aTrack.GetUserInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ParticleDefinition* ExExChProcessChanneling::
GetParticleDefinition(const G4Track& aTrack){
    return const_cast<G4ParticleDefinition*>(aTrack.GetParticleDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::
ComputeCrystalCharacteristic(const G4Track& aTrack){
    fIntegratedDensity->SetXPhysicalLattice(GetXPL(aTrack));
    fIntegratedDensity->InitializeTables();
    
    if(fFileCharacteristicsName != ""){
        G4String filename;
        fElectronDensity->InitializePhysicalLattice(GetXPL(aTrack));
        fNucleiDensity->InitializePhysicalLattice(GetXPL(aTrack));
        fPotentialEnergy->InitializePhysicalLattice(GetXPL(aTrack));
        fElectricField->InitializePhysicalLattice(GetXPL(aTrack));

        filename = fFileCharacteristicsName + "_pot.txt";
        G4cout << filename << G4endl;
        fPotentialEnergy->ReadFromECHARM(filename,CLHEP::eV);
        fPotentialEnergy->PrintOnFile("temppot.dat",GetXPL(aTrack));

        filename = fFileCharacteristicsName + "_efx.txt";
        G4cout << filename << G4endl;
        fElectricField->ReadFromECHARM(filename,CLHEP::eV/CLHEP::m);
        fElectricField->PrintOnFile("tempefx.dat",GetXPL(aTrack));

        filename = fFileCharacteristicsName + "_atd.txt";
        G4cout << filename << G4endl;
        fNucleiDensity->ReadFromECHARM(filename);

        filename = fFileCharacteristicsName + "_eld.txt";
        G4cout << filename << G4endl;
        fElectronDensity->ReadFromECHARM(filename);
        fIntegratedDensity->ReadFromFiles(fFileCharacteristicsName);
    }
    else{
        fPotentialEnergy->InitializePhysicalLattice(
                                                GetXPL(aTrack));
        fElectricField->InitializePhysicalLattice(GetXPL(aTrack));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
