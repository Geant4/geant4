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

void ExExChProcessChanneling::UpdateParameters(const G4Track& aTrack){
    
    if(fIntegratedDensity->HasBeenInitialized(GetXPhysicalLattice(aTrack))
       == false){
        ComputeCrystalCharacteristic(aTrack);
        G4cout << "ChannelingProcess::UpdatePositionMomentumDensity::";
        G4cout<<"fIntegratedDensity->Initialized" << G4endl;
    }
    
    UpdatePosition(aTrack);
    UpdateMomentum(aTrack);
    UpdateDensity(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::UpdatePosition(const G4Track& aTrack){
    
    if(GetInfo(aTrack)->GetPositionChanneledInitial().x() == DBL_MAX ||
       HasLatticeOnBoundaryPost(aTrack)){
        // when the particle enters the crystal the position in the channel
        //is randomly generated using a uniform distribution
        G4double vXposition = G4UniformRand() *
            GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
        
        //initial position in the channel is stored
        GetInfo(aTrack)->SetPositionChanneledInitial(G4ThreeVector(vXposition,
                                                                   0.,
                                                                   0.));
        
        //initial position for the measurement of apparent centrifugal force
        //is zero at crystal entrance
        GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(0.,0.,0.));
    }
    else{
        G4double vPositionX = GetInfo(aTrack)->GetPositionChanneled().x();
        
        //if the particle has been under channeling the position
        //for the measurement of the apparent centrifugal force is reset
        if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
            GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(0.,0.,0.));
        }
        else{
            //if the particle has not been under channeling the distance
            //between the new and the old point is computed for the evaluation
            //of the centrifugal potential acting on the particle
            vPositionX += (ComputePositionInTheCrystal(
                aTrack.GetStep()->GetPostStepPoint(),
                aTrack).x() - ComputePositionInTheCrystal(
                            aTrack.GetStep()->GetPreStepPoint(),
                            aTrack).x());
            GetInfo(aTrack)->SetPositionChanneled(G4ThreeVector(vPositionX,
                                                                0.,
                                                                0.));
        }
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::UpdateMomentum(const G4Track& aTrack){
    if(GetInfo(aTrack)->GetMomentumChanneledInitial().x() == DBL_MAX){
        // the first time it enter the crystal we take the momentum
        // for the post step which is the only one in the crystal
        G4ThreeVector vMomentum =
            ComputeMomentum(aTrack,aTrack.GetStep()->GetPostStepPoint());
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
        
        GetInfo(aTrack)->SetMomentumChanneledInitial(
                            GetInfo(aTrack)->GetMomentumChanneled());
    }
    else{
        // we take the PREVIOUS step point to compare,
        // otherwise the momentum is not computed correctly
        G4ThreeVector vMomentum =
            G4ThreeVector(GetInfo(aTrack)->GetMomentumChanneled().x(),
                          GetInfo(aTrack)->GetMomentumChanneled().y(),0.);
        
        vMomentum+=ComputeMomentum(aTrack,aTrack.GetStep()->GetPreStepPoint());
        
        GetInfo(aTrack)->SetMomentumChanneled(vMomentum);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChProcessChanneling::UpdateDensity(const G4Track& aTrack){
    
    
    if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
        
        G4double vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
        
        if(GetXPhysicalLattice(aTrack)->IsBent()){
            if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
                if(ParticleIsNegative(aTrack)){
                    vTransverseEnergy -=
                    ComputeCentrifugalEnergyMaximumVariation(aTrack).x();
                }
                else{
                    vTransverseEnergy +=
                    ComputeCentrifugalEnergyMaximumVariation(aTrack).x();
                }
            }
        }
        
        G4double vCharge = GetParticleDefinition(aTrack)->GetPDGCharge();
        G4double vNucleiDensity =
            fIntegratedDensity->GetIntegratedDensityNuclei(vTransverseEnergy,
                                            GetXPhysicalLattice(aTrack),
                                            G4int(vCharge));
        G4double vElectronDensity =
            fIntegratedDensity->GetIntegratedDensityElectron(vTransverseEnergy,
                                            GetXPhysicalLattice(aTrack),
                                            G4int(vCharge));
        
        G4double vLowerBoundNegative = 1.;
        G4double vLowerBoundPositive = 0.01;
        
        if(ParticleIsNegative(aTrack)){
            if(vNucleiDensity < vLowerBoundNegative)
            {vNucleiDensity = vLowerBoundNegative;}
            if(vElectronDensity < vLowerBoundNegative)
            {vElectronDensity = vLowerBoundNegative;}
        }
        else{
            if(vNucleiDensity < vLowerBoundPositive)
            {vNucleiDensity = vLowerBoundPositive;}
            if(vElectronDensity < vLowerBoundPositive)
            {vElectronDensity = vLowerBoundPositive;}
        }
        
        GetInfo(aTrack)->SetNucleiDensity(vNucleiDensity);
        GetInfo(aTrack)->SetElectronDensity(vElectronDensity);
    }
    else{
        ResetDensity(aTrack);
    }
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
ComputeChannelingOutgoingMomentum(const G4Track& aTrack){
    
    G4StepPoint* vStepPre = aTrack.GetStep()->GetPreStepPoint();
    G4StepPoint* vStepPost = aTrack.GetStep()->GetPostStepPoint();
    
    G4double vTotalEnergy = vStepPre->GetTotalEnergy();
    
    G4double vTransverseEnergyX = std::fabs(ComputeTransverseEnergy(aTrack).x());
    G4double vTransverseEnergyY = std::fabs(ComputeTransverseEnergy(aTrack).y());
    double vPotentialEnergyX = 0.;
    double vPotentialEnergyY = 0.;
    
    bool bExit = false;
    do{
        G4double vXposition = G4UniformRand() *
        GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
        
        GetInfo(aTrack)->SetPositionChanneledInitial(G4ThreeVector(vXposition,
                                                                   0.,
                                                                   0.));
        
        vPotentialEnergyX = ComputePotentialEnergy(aTrack).x();
        vPotentialEnergyY = ComputePotentialEnergy(aTrack).y();
        if(vPotentialEnergyX<=vTransverseEnergyX &&
           vPotentialEnergyY<=vTransverseEnergyY){
            bExit = true;
        }
    } while(bExit == false);

    vTransverseEnergyX-=vPotentialEnergyX;
    vTransverseEnergyY-=vPotentialEnergyY;

    G4double vChAngleX = std::pow(+ 2. * std::fabs(vTransverseEnergyX)
                             / vTotalEnergy , 0.5);
    G4double vChAngleY = std::pow(+ 2. * std::fabs(vTransverseEnergyY)
                             / vTotalEnergy , 0.5);
    
    G4double vPhi = 2. * ( G4UniformRand() - 0.5) * vChAngleX;
    G4double vTheta = 2. * ( G4UniformRand() - 0.5) * vChAngleY;
    
    G4ThreeVector vNewMomentum =
        G4ThreeVector(0.,0.,1.).rotate(G4ThreeVector(0,1,0),- vPhi)
        .rotate(G4ThreeVector(1,0,0),- vTheta);
    G4ThreeVector vPosition = ComputePositionInTheCrystal(vStepPost,aTrack);
    
    return GetXPhysicalLattice(aTrack)->
        ProjectMomentumVectorFromLatticeToWorld(vNewMomentum,vPosition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector ExExChProcessChanneling::
ComputeVolumeReflectionOutgoingMomentum(const G4Track& aTrack){
    
    G4StepPoint* vStep = aTrack.GetStep()->GetPostStepPoint();
    
    G4double vVrAngle = 0.;
    
    if(GetXPhysicalLattice(aTrack)->IsBent()) {
        G4double vRadiusX = GetXPhysicalLattice(aTrack)->
            GetCurvatureRadius().x();
        
        G4double vTotalEnergy = vStep->GetTotalEnergy();
        
        G4double vEnergyMax =
            std::fabs(ComputeCriticalEnergyMaximum(aTrack)
                - ComputeCriticalEnergyMinimum(aTrack));
        
        G4double vEnergyRMS =
            std::fabs(ComputeCentrifugalEnergyMaximumVariation(aTrack).x());
        
        G4double vTransverseEnergy =
            vEnergyMax + (G4UniformRand() * std::fabs(vEnergyRMS) );
        
        vVrAngle = - std::fabs(vRadiusX)/vRadiusX *
            std::pow(+ 2. * std::fabs(vTransverseEnergy) / vTotalEnergy , 0.5);
        
        if(ParticleIsNegative(aTrack)){
            vVrAngle *= 0.8; // = see PLB 681 (2009) 233
        }
        else{
            vVrAngle *= 1.4;
        }
        
        G4ThreeVector vMomentumChanneled =
            GetInfo(aTrack)->GetMomentumChanneled();
        G4double vAngleRatio =
            (vMomentumChanneled.x()/vTotalEnergy)/ComputeCriticalAngle(aTrack);
        
        if(std::fabs(vAngleRatio)<1.5){
            vVrAngle *= (-(std::fabs(vAngleRatio) - 1.5)/3.);
        }
    }
    
    G4double vOmega = GetXPhysicalLattice(aTrack)->GetLatticeAngles().y();
    G4double vPhi = vVrAngle * std::cos(vOmega);
    G4double vTheta = vVrAngle * std::sin(vOmega);
    
    G4ThreeVector vNewMomentum =
        aTrack.GetMomentum().unit()
        .rotate(G4ThreeVector(0.,1.,0.), - vPhi)
        .rotate(G4ThreeVector(1.,0.,0.), -vTheta);
    
    return vNewMomentum.unit();
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
    
    
    if(GetXPhysicalLattice(aTrack)->IsBent() == false){
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

G4double ExExChProcessChanneling::
ComputeDistanceWhereParticleTangentToBentPlane(const G4Track& aTrack){
    G4ThreeVector vPositionPre =
        ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
    G4ThreeVector vMomentumPre =
        GetXPhysicalLattice(aTrack)->
            ProjectMomentumVectorFromWorldToLattice(
                        aTrack.GetStep()->GetPreStepPoint()->GetMomentum(),
                        vPositionPre);
    
    G4double vDeltaProportion = 1.;
    
    if((vMomentumPre.x())!=0.){
        vDeltaProportion = std::fabs(vMomentumPre.unit().x());
    }
    G4double vDeltaPosition = vDeltaProportion*
        GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x();
    
    return std::abs(vDeltaPosition);
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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool ExExChProcessChanneling::
IsUnderCoherentEffect(const G4Track& aTrack){
    //----------------------------------------
    // check if the particle momentum
    // transverse to the (h,k,l) plane
    // is small enough to permit channeling
    //----------------------------------------
    
    UpdateParameters(aTrack);
    
    G4double vEnergyMax = ComputeCriticalEnergyMaximum(aTrack);
    
    G4double vTransverseEnergy = ComputeTransverseEnergy(aTrack).x();
    
    if(GetXPhysicalLattice(aTrack)->IsBent() == false){
        if(vTransverseEnergy <= vEnergyMax){
            GetInfo(aTrack)->SetCoherentEffect(1);
            // the particle is in channeling
            return true;
        }
    }
    else{
        G4ThreeVector vPositionInTheCrystal =
        GetInfo(aTrack)->GetPositionChanneled()
            + GetInfo(aTrack)->GetPositionChanneledInitial();
        vTransverseEnergy += std::fabs(ComputeCentrifugalEnergy(aTrack,
                                                vPositionInTheCrystal).x());
        if(vTransverseEnergy <= vEnergyMax){
            // the particle is in channeling
            GetInfo(aTrack)->SetCoherentEffect(1);
            return true;
        }
        else{
            G4bool bNotBoundary = ParticleIsNotOnBoundary(aTrack);
            G4bool bTangentToPlane = ParticleIsTangentToBentPlane(aTrack);

            if(bTangentToPlane == true &&
               bNotBoundary == true &&
               GetInfo(aTrack)->HasBeenUnderCoherentEffect() != 2){
                    // the particle is in volume reflection
                    GetInfo(aTrack)->SetCoherentEffect(2);
                    return true;
                }
            }
    }
    
    // the particle is not under coherent effect
    GetInfo(aTrack)->SetCoherentEffect(0);
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double ExExChProcessChanneling::
GetChannelingMeanFreePath(const G4Track& aTrack){
    //----------------------------------------
    // return the channeling MFP
    //----------------------------------------
    
    G4double vMFPosc = ComputeOscillationPeriod(aTrack);
    
    if(GetInfo(aTrack)->GetNucleiDensity() < 1.){
        vMFPosc /= GetInfo(aTrack)->GetNucleiDensity();
    }
    
    G4double vMFP = vMFPosc * 2.;
    
    if(GetXPhysicalLattice(aTrack)->IsBent()){
        G4double vMFPVR =
            ComputeDistanceWhereParticleTangentToBentPlane(aTrack);
        
        if((std::fabs(vMFPVR) < vMFP) && (std::fabs(vMFPVR) > (0.5 * vMFPosc))){
            vMFP = vMFPVR;
        }
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
        return GetChannelingMeanFreePath(aTrack);
    }
    else{
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
    
    G4bool bIsUnderCoherentEffect = false;
    
    if((HasLattice(aTrack) == true) &&
       (ParticleIsNotOnBoundaryPost(aTrack) == true)){
        bIsUnderCoherentEffect = IsUnderCoherentEffect(aTrack);
        if(bIsUnderCoherentEffect){
            if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
                // if the particle is in channeling it gives the direction
                // of the lattice to the particle momentum
                G4ThreeVector vPosition =
                    ComputePositionInTheCrystal(
                            aTrack.GetStep()->GetPostStepPoint(),aTrack);
                aParticleChange.ProposeMomentumDirection(
                    GetXPhysicalLattice(aTrack)->
                    GetLatticeDirection(vPosition).unit());
            }
            else if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 2){
                // if the particle is in VR it gives a kick
                // to the opposite site of the bending to the particle
                aParticleChange.ProposeMomentumDirection(
                            ComputeVolumeReflectionOutgoingMomentum(aTrack));
                GetInfo(aTrack)->SetCoherentEffect(0);
                ResetDensity(aTrack);
            }
        }
    }
    else{
        // if the volume has no lattice it resets the density factors
        ResetDensity(aTrack);
    }
    
    if( (bIsUnderCoherentEffect == false && (HasLattice(aTrack) == true) )
       || (HasLatticeOnBoundaryPre(aTrack) == true) ) {
        // if has been under coherent effect but now it is not,
        // the outgoing momentum is evaluated starting from the current position
        if(GetInfo(aTrack)->HasBeenUnderCoherentEffect() == 1){
            aParticleChange.ProposeMomentumDirection(
                            ComputeChannelingOutgoingMomentum(aTrack));
        }
        
        // If is not under coherent effect sets coherent effect to zero
        // and resets the density factors after the outgoing angle
        // has been evaluated
        GetInfo(aTrack)->SetCoherentEffect(0);
        ResetDensity(aTrack);
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
    //while(!getchar());
    return vTransverseEnergy;
}

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
        fPotentialEnergy->GetEC(GetInfo(aTrack)->GetPositionChanneledInitial(),
                                GetXPhysicalLattice(aTrack));
    
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
    
    G4ThreeVector vEnergyVariation =
        G4ThreeVector(vTotalEnergy * vPositionX /
                      GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x(),
                      0.,
                      0.);
    
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
    
    G4ThreeVector vMomentum =
        GetXPhysicalLattice(aTrack)->
            ProjectMomentumVectorFromWorldToLattice(aTrack.GetMomentum(),
                                                    vPosition);
    
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
    
    G4ThreeVector vEnergyVariation = G4ThreeVector(vTotalEnergy *
                    GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod() /
                    GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x(),
                    0.,
                    0.);
    
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
            - fPotentialEnergy->GetMinimum(GetXPhysicalLattice(aTrack));
    }
    else{
        vCriticalEnergy =
            + fPotentialEnergy->GetMaximum(GetXPhysicalLattice(aTrack));
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
            - fPotentialEnergy->GetMaximum(GetXPhysicalLattice(aTrack));
    }
    else{
        vCriticalEnergy =
            + fPotentialEnergy->GetMinimum(GetXPhysicalLattice(aTrack));
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
        GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
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
    G4double vCriticalRadius =
        vTotalEnergy / fElectricField->GetMaximum(GetXPhysicalLattice(aTrack));
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
        0.5 * GetXPhysicalLattice(aTrack)->ComputeInterplanarPeriod();
    
    G4double vCentreX = vInterplanarPeriodHalf;
    
    if(GetXPhysicalLattice(aTrack)->IsBent()){
        G4double vTotalEnergy =
            aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
        
        G4double vPotentialWellDepth =
            ComputeCriticalEnergyMaximum(aTrack)
            - (ComputeCriticalEnergyMinimum(aTrack));
        
        vCentreX *= (1. - 0.5 * vTotalEnergy /
                     vPotentialWellDepth /
                     GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x() *
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
GetXPhysicalLattice(const G4Track& aTrack){
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
       aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
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
       aTrack.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
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

G4bool ExExChProcessChanneling::
ParticleIsTangentToBentPlane(const G4Track& aTrack){
    G4ThreeVector vPositionPre =
        ComputePositionInTheCrystal(aTrack.GetStep()->GetPreStepPoint(),aTrack);
    G4ThreeVector vMomentumPre =
        GetXPhysicalLattice(aTrack)->
        ProjectMomentumVectorFromWorldToLattice(aTrack.GetStep()->
                GetPreStepPoint()->GetMomentum(),vPositionPre);
    
    G4ThreeVector vPositionPost =
        ComputePositionInTheCrystal(aTrack.GetStep()->
                GetPostStepPoint(),aTrack);
    G4ThreeVector vMomentumPost = GetXPhysicalLattice(aTrack)->
        ProjectMomentumVectorFromWorldToLattice(aTrack.GetStep()->
                GetPostStepPoint()->GetMomentum(),vPositionPost);
    
    if(vMomentumPost.x()<0. &&
       vMomentumPre.x()>0. &&
       GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x() < 0.){
        return true;
    }
    if(vMomentumPost.x()>0. &&
       vMomentumPre.x()<0. &&
       GetXPhysicalLattice(aTrack)->GetCurvatureRadius().x() > 0.){
        return true;
    }
    
    return false;
}

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
    fIntegratedDensity->SetXPhysicalLattice(GetXPhysicalLattice(aTrack));
    fIntegratedDensity->InitializeTables();
    
    if(fFileCharacteristicsName != ""){
        G4String filename;
        
        fPotentialEnergy->ReadFromECHARM(filename =
                                         fFileCharacteristicsName + "_pot.txt");
        //fElectricField->ReadFromECHARM("efx.txt");
        fElectricField->InitializePhysicalLattice(GetXPhysicalLattice(aTrack));
        fIntegratedDensity->ReadFromFiles(fFileCharacteristicsName);
    }
    else{
        fPotentialEnergy->InitializePhysicalLattice(
                                                GetXPhysicalLattice(aTrack));
        fElectricField->InitializePhysicalLattice(GetXPhysicalLattice(aTrack));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
