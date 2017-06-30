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

#include "G4Channeling.hh"

#include "Randomize.hh"

#include "G4ChannelingTrackData.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include "G4LambdacPlus.hh"

G4Channeling::G4Channeling():
G4VDiscreteProcess("channeling"),
fChannelingID(-1),
fTimeStepMin(0.),
fTimeStepMax(0.),
fTransverseVariationMax(2.E-2 * CLHEP::angstrom),
k010(G4ThreeVector(0.,1.,0.)){
    fChannelingID = G4PhysicsModelCatalog::GetIndex("channeling");
    if(fChannelingID == -1){
        fChannelingID = G4PhysicsModelCatalog::Register("channeling");
    }
    fSpin = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Channeling::~G4Channeling(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingTrackData* G4Channeling::GetTrackData(const G4Track& aTrack){
    G4ChannelingTrackData* trackdata =
        (G4ChannelingTrackData*)(aTrack.GetAuxiliaryTrackInformation(fChannelingID));
    if(trackdata == nullptr){
        trackdata = new G4ChannelingTrackData();
        aTrack.SetAuxiliaryTrackInformation(fChannelingID,trackdata);
    }
    return trackdata;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::GetEF(const G4Track& aTrack,
                         G4ThreeVector& pos,
                         G4ThreeVector& out){
    out = G4ThreeVector((GetMatData(aTrack)->GetEFX()->GetEC(pos)),
                        (GetMatData(aTrack)->GetEFY()->GetEC(pos)),
                        0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::PosToLattice(G4StepPoint* step,G4ThreeVector& pos){
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(step->GetTouchable());

    pos -= theTouchable->GetTranslation();
    pos = ((*theTouchable->GetRotation()).inverse())(pos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4Channeling::UpdateParameters(const G4Track& aTrack){

    G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
    
    G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
    G4StepPoint* preStepPoint = aTrack.GetStep()->GetPreStepPoint();

    G4ThreeVector posPost = postStepPoint->GetPosition();
    aLCV->RotateToLattice(posPost);
    G4ThreeVector posPre = preStepPoint->GetPosition();
    aLCV->RotateToLattice(posPre);
    
    G4double integrationLimit = fabs(posPost.z() - posPre.z());
    
    if(integrationLimit > 0.){
        //----------------------------------------
        // Check if the crystal is bent
        //----------------------------------------
        G4bool isBent = GetMatData(aTrack)->IsBent();

        //----------------------------------------
        // Get the momentum in the world reference
        // frame and rotate to the solid reference frame
        //----------------------------------------
        G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
        G4ThreeVector momWorld = aTrack.GetStep()->GetPreStepPoint()->GetMomentum();
        G4ThreeVector mom = (*theTouchable->GetRotation())(momWorld);
        
        //----------------------------------------
        // Get the momentum in the solid reference
        // frame and rotate to the crystal reference frame
        //----------------------------------------
        aLCV->RotateToLattice(mom);

        //----------------------------------------
        // Get the momentum in the crystal reference
        // frame and rotate to the reference frame
        // solidal to the bent planes
        //----------------------------------------
        if(isBent){
            PosToLattice(preStepPoint,posPre);
            G4ThreeVector axis010 = (*theTouchable->GetRotation())(k010);
            mom.rotate(axis010,-posPre.z()/GetMatData(aTrack)->GetBR(posPre).x());
        }

        //----------------------------------------
        // Take the position stored in the track data.
        // If the particle enters the crystal,
        // the position in the channel is randomly
        // generated using a uniform distribution
        //----------------------------------------
        G4ThreeVector pos;
        if(GetTrackData(aTrack)->GetPosCh().x() == DBL_MAX){
            G4double posX = G4UniformRand() * GetMatData(aTrack)->GetPot()->GetIntSp(0);
            G4double posY = G4UniformRand() * GetMatData(aTrack)->GetPot()->GetIntSp(1);
            pos = G4ThreeVector(posX,posY,0.);
        }
        else{
            pos = GetTrackData(aTrack)->GetPosCh();
        }

        G4double step=0., stepTot=0.;
        G4double nud =0., eld    =0.;
        G4double efx =0., efy    =0.;
        G4double nud_temp =0., eld_temp    =0.;

        G4double beta = aTrack.GetVelocity()/CLHEP::c_light;
        G4double Z = GetParticleDefinition(aTrack)->GetPDGCharge();
        
        const G4double oneSixth = 1./6.;
        G4ThreeVector posk1,posk2,posk3,posk4,posk5,posk6;
        G4ThreeVector momk1,momk2,momk3,momk4,momk5,momk6;
        G4ThreeVector pos_temp, efxy;

        do{
            //----------------------------------------
            // Limit the variable step length for the
            // integration via the selected algorithm
            // and update variables for the integration
            //----------------------------------------

            UpdateIntegrationStep(aTrack,mom,step);
            if(step + stepTot > integrationLimit){
                step = integrationLimit - stepTot;
            }

            //----------------------------------------
            // Function integration algorithm
            // 4th Order Runge-Kutta
            //----------------------------------------
            
            GetEF(aTrack,pos,efxy);
            posk1 = step / mom.z() * mom;
            momk1 = step / beta * Z * efxy;
            if(isBent) momk1.setX(momk1.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos)).x());
            
            GetEF(aTrack,pos_temp = pos + posk1 * 0.5,efxy);
            posk2 = step / mom.z() * (mom + momk1 * 0.5);
            momk2 = step / beta * Z * efxy;
            if(isBent) momk2.setX(momk2.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos_temp)).x());

            GetEF(aTrack,pos_temp = pos + posk2 * 0.5,efxy);
            posk3 = step / mom.z() * (mom + momk2 * 0.5);
            momk3 = step / beta * Z * efxy;
            if(isBent) momk3.setX(momk3.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos_temp)).x());
            
            GetEF(aTrack,pos_temp = pos + posk3,efxy);
            posk4 = step / mom.z() * (mom + momk3);
            momk4 = step / beta * Z * efxy;
            if(isBent) momk4.setX(momk4.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos_temp)).x());

            pos = pos + oneSixth * (posk1 + 2.*posk2 + 2.*posk3 + posk4);
            mom = mom + oneSixth * (momk1 + 2.*momk2 + 2.*momk3 + momk4);
            
            //----------------------------------------
            // Function integration algorithm
            // 2th Order Velocity-Verlet
            //----------------------------------------
           
            /*
            GetEF(aTrack,pos,efxy);
            posk1 = pos + (step * 0.5 / mom.z()) * mom;
            //momk1 = mom + step * 0.5 / betaZ * efxy;
            momk1 = mom;
            if(isBent) momk1.setX(momk1.x() - step * 0.5 * mom.z() * beta / (GetMatData(aTrack)->GetBR(pos)).x());
            
            GetEF(aTrack,posk1,efxy);
            pos = pos + (step / momk1.z()) * momk1;
            //mom = mom + step / betaZ * efxy;
            mom = mom;
            if(isBent) mom.setX(mom.x() - step * mom.z() * beta / (GetMatData(aTrack)->GetBR(posk1)).x());
            */

            //----------------------------------------
            // Update the total step and the electron
            // and nuclei density experienced by
            // the particle during its motion
            //----------------------------------------

            stepTot += step;

            nud_temp = GetMatData(aTrack)->GetNuD()->GetEC(pos);
            eld_temp = GetMatData(aTrack)->GetElD()->GetEC(pos);
            
            if(nud_temp < 0.) {nud_temp = 0.;}
            if(eld_temp < 0.) {eld_temp = 0.;}

            nud += (step * nud_temp);
            eld += (step * eld_temp);

            efx += (step * GetMatData(aTrack)->GetEFX()->GetEC(pos));
            efy += (step * GetMatData(aTrack)->GetEFY()->GetEC(pos));

        } while(stepTot<integrationLimit);
        
        nud /= stepTot;
        eld /= stepTot;

        if(nud < 1.E-10) {nud = 1.E-10;}
        if(eld < 1.E-10) {eld = 1.E-10;}
        
        GetTrackData(aTrack)->SetNuD(nud);
        GetTrackData(aTrack)->SetElD(eld);

        GetTrackData(aTrack)->SetEFX(efx);
        GetTrackData(aTrack)->SetEFY(efy);
        
        GetTrackData(aTrack)->SetMomCh(mom);
        GetTrackData(aTrack)->SetPosCh(pos);
        return true;
    }
    else{
        return false;
    }
    
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4Channeling::
UpdateIntegrationStep(const G4Track& aTrack,
                      G4ThreeVector& mom,
                      G4double& step){
    
    if(mom.x() != 0.0 || mom.y() != 0.0){
        double xy2 = mom.x() * mom.x() + mom.y()*mom.y();
        
        if(xy2!=0.){
            step = std::fabs(fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy() / std::pow(xy2,0.5));
            if(step < fTimeStepMin) step = fTimeStepMin;
            else{
                fTimeStepMax = sqrt( fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy()
                                    / fabs(GetMatData(aTrack)->GetEFX()->GetMax()));
                
                if(step > fTimeStepMax) step = fTimeStepMax;
            }
        }
        else{
            step = fTimeStepMin;
        }
        
        return true;
    }
    else{
        step = fTimeStepMin;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Channeling::
GetMeanFreePath(const G4Track& aTrack,
                G4double, // previousStepSize
                G4ForceCondition* condition){
    
    //----------------------------------------
    // the condition is forced to check if
    // the volume has a lattice at each step.
    // if it hasn't, return DBL_MAX
    //----------------------------------------
    
    *condition = Forced;
    
    G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
    G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();
    
    if(G4LogicalCrystalVolume::IsLattice(aLV) == true &&
       G4LogicalCrystalVolume::IsLattice(aNLV) == true){
        G4double osc_per = GetOscillationPeriod(aTrack);
        fTimeStepMin = osc_per * 2.E-4;
        return osc_per * 0.01;
    }
    else{
        GetTrackData(aTrack)->Reset();
        return DBL_MAX;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4Channeling::
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
    G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
    G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();
    
    
    if(G4LogicalCrystalVolume::IsLattice(aLV) == true &&
       G4LogicalCrystalVolume::IsLattice(aNLV) == true){

        G4bool bModifiedTraj = UpdateParameters(aTrack);

        if(bModifiedTraj==true){
            //----------------------------------------
            // Get the momentum in the reference frame
            // solidal to the bent planes and rotate
            // to the reference frame
            //----------------------------------------
            G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
            G4ThreeVector momCh = GetTrackData(aTrack)->GetMomCh();

            G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
            G4TouchableHistory* theTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());

            if(GetMatData(aTrack)->IsBent()){
                G4ThreeVector posPost = postStepPoint->GetPosition();
                PosToLattice(postStepPoint,posPost);
                G4ThreeVector axis010 = (*theTouchable->GetRotation())(k010);
                momCh.rotate(axis010,posPost.z()/GetMatData(aTrack)->GetBR(posPost).x());
            }
            
            //----------------------------------------
            // Get the momentum in the crystal reference
            // frame and rotate to the solid reference frame
            //----------------------------------------

            aLCV->RotateToSolid(momCh);
            
            //----------------------------------------
            // Get the momentum in the solid reference
            // frame and rotate to the world reference frame
            //----------------------------------------
            G4ThreeVector mom = ((*theTouchable->GetRotation()).inverse())(momCh);

            aParticleChange.ProposeMomentumDirection(mom.unit());
            aParticleChange.ProposePolarization(fSpin);
        }
    }
    else{
        // if the volume has no lattice it resets the density factors
        GetTrackData(aTrack)->Reset();
    }
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
