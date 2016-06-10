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

#include "XWrapperContinuousDiscreteProcess.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"

#include "G4SystemOfUnits.hh"


XWrapperContinuousDiscreteProcess::
XWrapperContinuousDiscreteProcess(const G4String& aName)
:G4VContinuousDiscreteProcess(aName){
    if (verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    bNucleiOrElectronFlag = +0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XWrapperContinuousDiscreteProcess::
XWrapperContinuousDiscreteProcess(const G4String& aName,
                                  G4VContinuousDiscreteProcess* toRegister)
:G4VContinuousDiscreteProcess(aName){
    fRegisteredProcess = toRegister;
    if (verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    bNucleiOrElectronFlag = +0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XWrapperContinuousDiscreteProcess::
XWrapperContinuousDiscreteProcess(const G4String& aName,
                                  G4ProcessType aProcessType)
:G4VContinuousDiscreteProcess(aName,aProcessType){
    if (verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    bNucleiOrElectronFlag = +0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XWrapperContinuousDiscreteProcess::~XWrapperContinuousDiscreteProcess(){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XWrapperContinuousDiscreteProcess::
XWrapperContinuousDiscreteProcess(XWrapperContinuousDiscreteProcess& right):
G4VContinuousDiscreteProcess(right){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XWrapperContinuousDiscreteProcess::
RegisterProcess(G4VContinuousDiscreteProcess* toRegister){
    fRegisteredProcess = toRegister;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XWrapperContinuousDiscreteProcess::
RegisterProcess(G4VContinuousDiscreteProcess* toRegister,
                G4int flag){
    fRegisteredProcess = toRegister;
    bNucleiOrElectronFlag = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XWrapperContinuousDiscreteProcess::SetNucleiOrElectronFlag(G4int flag){
    bNucleiOrElectronFlag = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int XWrapperContinuousDiscreteProcess::GetNucleiOrElectronFlag(){
    return bNucleiOrElectronFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XWrapperContinuousDiscreteProcess::GetDensity(const G4Track& aTrack){
    //Retrieve nuclei and electron density
    //from ExExChParticleUserInfo object
    ExExChParticleUserInfo* chanInfo =
        (ExExChParticleUserInfo*) aTrack.GetUserInformation();
    
    G4double vDensity = 1.;
    
    if(chanInfo){
        if(bNucleiOrElectronFlag == +1){
            vDensity = chanInfo->GetNucleiDensity();
        }
        else if(bNucleiOrElectronFlag == -1){
            vDensity = chanInfo->GetElectronDensity();
        }
        else{
            vDensity = (chanInfo->GetNucleiDensity()
                        + chanInfo->GetElectronDensity())/2.;
        }
    }
    else {
        G4cout << G4endl << "XWrapperContinuousDiscreteProcess::";
        G4cout << "ERROR - no ExExChParticleUserInfo object Detected";
        G4cout << G4endl;
    }
    
    return vDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XWrapperContinuousDiscreteProcess::
GetDensityPreviousStep(const G4Track& aTrack){
    //Retrieve nuclei and electron density
    //from ExExChParticleUserInfo object
    ExExChParticleUserInfo* chanInfo =
        (ExExChParticleUserInfo*) aTrack.GetUserInformation();
    
    G4double vDensityPreviousStep = 1.;
    
    if(chanInfo){
        if(bNucleiOrElectronFlag == +1){
            vDensityPreviousStep = chanInfo->GetNucleiDensityPreviousStep();
        }
        else if(bNucleiOrElectronFlag == -1){
            vDensityPreviousStep = chanInfo->GetElectronDensityPreviousStep();
        }
        else{
            vDensityPreviousStep =
                (chanInfo->GetNucleiDensityPreviousStep()
                 + chanInfo->GetElectronDensityPreviousStep())/2.;
        }
    }
    else {
        G4cout << G4endl << "XWrapperContinuousDiscreteProcess::";
        G4cout << "ERROR - no ExExChParticleUserInfo";
        G4cout << "object Detected" << G4endl;
    }
    
    return vDensityPreviousStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XWrapperContinuousDiscreteProcess::StartTracking(G4Track* aTrack){
    fRegisteredProcess->StartTracking(aTrack);
    currentInteractionLength = -1.0;
    theNumberOfInteractionLengthLeft = -1.0;
    theInitialNumberOfInteractionLength = -1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XWrapperContinuousDiscreteProcess::
GetMeanFreePath(const G4Track&,
                G4double, //previousStepSize,
                G4ForceCondition*){
    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XWrapperContinuousDiscreteProcess::
PostStepGetPhysicalInteractionLength (const G4Track &aTrack,
                                      G4double previousStepSize,
                                      G4ForceCondition *condition){
    
    G4double vDensity = GetDensity(aTrack);
    G4double vDensityPreviousStep = GetDensityPreviousStep(aTrack);
    
    if ( (previousStepSize < 0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
        // beginning of tracking (or just after DoIt of this process)
        ResetNumberOfInteractionLengthLeft();
    } else if ( previousStepSize > 0.0) {
        // subtract NumberOfInteractionLengthLeft
        SubtractNumberOfInteractionLengthLeft(previousStepSize
                                              * vDensityPreviousStep);
    } else {
        // zero step DO NOTHING
    }
    
    G4double regIntLength =
        fRegisteredProcess->PostStepGetPhysicalInteractionLength(aTrack,
                                previousStepSize * vDensityPreviousStep,
                                condition);
    G4double regIntNumber =
        fRegisteredProcess->GetNumberOfInteractionLengthLeft();
    currentInteractionLength = regIntLength / regIntNumber;
    theNumberOfInteractionLengthLeft = regIntNumber;
    
    currentInteractionLength =  theNumberOfInteractionLengthLeft
        * currentInteractionLength;
    currentInteractionLength /= vDensity;
    return currentInteractionLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XWrapperContinuousDiscreteProcess::
AlongStepGetPhysicalInteractionLength (const G4Track& aTrack,
                                       G4double  previousStepSize,
                                       G4double  currentMinimumStep,
                                       G4double& currentSafety,
                                       G4GPILSelection* selection){
    G4double vDensityPreviousStep = GetDensityPreviousStep(aTrack);

    return fRegisteredProcess->AlongStepGetPhysicalInteractionLength(aTrack,
                                    previousStepSize * vDensityPreviousStep,
                                    currentMinimumStep,
                                    currentSafety,
                                    selection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* XWrapperContinuousDiscreteProcess::
PostStepDoIt(const G4Track& aTrack,
             const G4Step& aStep ){
    G4double vDensity = GetDensity(aTrack);
    G4double vStepLengthSaved = aStep.GetStepLength();
    const_cast<G4Step&>(aStep).SetStepLength(aStep.GetStepLength() * vDensity);
    pParticleChange = fRegisteredProcess->PostStepDoIt(aTrack, aStep);
    const_cast<G4Step&>(aStep).SetStepLength(vStepLengthSaved);
    return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* XWrapperContinuousDiscreteProcess::
AlongStepDoIt(const G4Track& aTrack,
              const G4Step& aStep ){
    G4double vDensity = GetDensity(aTrack);
    G4double vStepLengthSaved = aStep.GetStepLength();
    const_cast<G4Step&>(aStep).SetStepLength(aStep.GetStepLength() * vDensity);
    pParticleChange = fRegisteredProcess->AlongStepDoIt(aTrack, aStep);
    const_cast<G4Step&>(aStep).SetStepLength(vStepLengthSaved);
    return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XWrapperContinuousDiscreteProcess::
GetContinuousStepLimit(const G4Track& ,
                       G4double,
                       G4double,
                       G4double& ){
    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
