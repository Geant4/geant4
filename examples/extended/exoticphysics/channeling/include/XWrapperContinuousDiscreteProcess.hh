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

#ifndef XWrapperContinuousDiscreteProcess_h
#define XWrapperContinuousDiscreteProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "ExExChParticleUserInfo.hh"
#include "G4ParticleChangeForNothing.hh"

class G4Material;

class XWrapperContinuousDiscreteProcess : public G4VContinuousDiscreteProcess
{
public:
    
    XWrapperContinuousDiscreteProcess(const G4String& processName
                                      ="XWrapperContinuousDiscreteProcess" );
    XWrapperContinuousDiscreteProcess(const G4String& processName,
                                      G4VContinuousDiscreteProcess*);
    XWrapperContinuousDiscreteProcess(const G4String& processName,
                                      G4ProcessType);
    
    G4int ItHasToWork(const G4Track&);
    
    virtual ~XWrapperContinuousDiscreteProcess();
    
public:
    void RegisterProcess(G4VContinuousDiscreteProcess*);
    void RegisterProcess(G4VContinuousDiscreteProcess*,G4int,G4int aBool = 0);

    G4VContinuousDiscreteProcess* GetRegisteredProcess()
    {return fRegisteredProcess;};

    G4double GetDensity(const G4Track&);
    G4double GetDensityPreviousStep(const G4Track&);
    
    void SetNucleiOrElectronFlag(G4int);
    G4int GetNucleiOrElectronFlag();
    
private:
    // hide assignment operator as private
    XWrapperContinuousDiscreteProcess(XWrapperContinuousDiscreteProcess&);
    XWrapperContinuousDiscreteProcess& operator=(
                            const XWrapperContinuousDiscreteProcess& right);
    
    //private data members
    G4int bBothOrCrystalOrDetectorPhysics;
    G4int bNucleiOrElectronFlag;
    //Decide whether to use nuclei (+1) or electron (-1)
    //or both (0) density to change parameters
    G4VContinuousDiscreteProcess* fRegisteredProcess;
    
    const G4Step theStepCopy;
    G4ParticleChangeForNothing* fParticleChangeForNothing;
    /////////////////////////////////////////////////////////
    /////////////////// GEANT4 PROCESS METHODS //////////////
    /////////////////////////////////////////////////////////
public:
    // DO IT
    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );
    virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step& );
    virtual G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
            const G4Step& aStep){
        return fRegisteredProcess->AtRestDoIt(aTrack,aStep);
    }

    // GPIL
    virtual G4double PostStepGetPhysicalInteractionLength (
                            const G4Track&,
                            G4double,
                            G4ForceCondition*);
    virtual G4double AlongStepGetPhysicalInteractionLength (
                            const G4Track&,
                            G4double,
                            G4double,
                            G4double&,
                            G4GPILSelection*);
    virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& aTrack,
            G4ForceCondition* condition){
        return fRegisteredProcess->AtRestGetPhysicalInteractionLength(aTrack,
               condition);
    };

    // GENERAL
    virtual void StartTracking(G4Track* aTrack);
    virtual void EndTracking() {fRegisteredProcess->EndTracking();};
    virtual void SetProcessManager(const G4ProcessManager* aPM)
    {fRegisteredProcess->SetProcessManager(aPM);};
    virtual const G4ProcessManager* GetProcessManager()
    {return fRegisteredProcess->GetProcessManager();};
    virtual void ResetNumberOfInteractionLengthLeft()
    {fRegisteredProcess->ResetNumberOfInteractionLengthLeft();};
    virtual void DumpInfo() const {fRegisteredProcess->DumpInfo();};
    
    
    virtual void SetMasterProcess(G4VProcess* masterP) {
        fRegisteredProcess->SetMasterProcess(
                static_cast<XWrapperContinuousDiscreteProcess*>(masterP)
                            ->fRegisteredProcess);
    };
    
    virtual void BuildWorkerPhysicsTable(const G4ParticleDefinition& aPD)
    {fRegisteredProcess->BuildWorkerPhysicsTable(aPD);};
    virtual void PrepareWorkerPhysicsTable(const G4ParticleDefinition& aPD)
    {fRegisteredProcess->PrepareWorkerPhysicsTable(aPD);};

    
    // PHYSICS TABLE
    virtual G4bool IsApplicable(const G4ParticleDefinition& aPD)
    {return fRegisteredProcess->IsApplicable(aPD);};
    virtual void BuildPhysicsTable(const G4ParticleDefinition& aPD)
    {fRegisteredProcess->BuildPhysicsTable(aPD);};
    virtual void PreparePhysicsTable(const G4ParticleDefinition& aPD)
    {fRegisteredProcess->PreparePhysicsTable(aPD);};
    virtual G4bool StorePhysicsTable(const G4ParticleDefinition* aPD,
                                     const G4String& aString,
                                     G4bool aBool)
    {return fRegisteredProcess->StorePhysicsTable(aPD,aString,aBool);};
    virtual G4bool RetrievePhysicsTable(const G4ParticleDefinition* aPD,
                                        const G4String& aString,
                                        G4bool aBool)
    {return fRegisteredProcess->RetrievePhysicsTable(aPD,aString,aBool);};
protected:
    // MFP
    virtual G4double GetMeanFreePath(const G4Track&,
                                     G4double,
                                     G4ForceCondition* );

protected:
    virtual G4double GetContinuousStepLimit(const G4Track&,
                                            G4double,
                                            G4double,
                                            G4double&);
    
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    
};

#endif
