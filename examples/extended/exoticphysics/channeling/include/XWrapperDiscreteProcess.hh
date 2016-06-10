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

#ifndef XWrapperDiscreteProcess_h
#define XWrapperDiscreteProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "ExExChParticleUserInfo.hh"
#include "G4ParticleChangeForNothing.hh"

class G4Material;

class XWrapperDiscreteProcess : public G4VDiscreteProcess
{
public:
    
    XWrapperDiscreteProcess(const G4String& processName =
                                        "XWrapperDiscreteProcess" );
    XWrapperDiscreteProcess(const G4String& , G4VDiscreteProcess*);
    
    G4int ItHasToWork(const G4Track&);

    virtual ~XWrapperDiscreteProcess();
    
public:
    void RegisterProcess(G4VDiscreteProcess*);
    void RegisterProcess(G4VDiscreteProcess*,G4int,G4int aBool = 0);
    
    G4double GetDensity(const G4Track&);
    G4double GetDensityPreviousStep(const G4Track&);
    
    void SetNucleiOrElectronFlag(G4int);
    G4int GetNucleiOrElectronFlag();
    
private:
    // hide assignment operator as private
    XWrapperDiscreteProcess(XWrapperDiscreteProcess&);
    XWrapperDiscreteProcess& operator=(const XWrapperDiscreteProcess&);
    
    //private data members
    ///Decide whether to use nuclei (+1) or electron (-1) or both (0)
    //density to change parameter
    G4int bBothOrCrystalOrDetectorPhysics;
    G4int bNucleiOrElectronFlag;
    G4VDiscreteProcess* fRegisteredProcess;
    G4ParticleChangeForNothing* fParticleChangeForNothing;

    
    /////////////////////////////////////////////////////////
    /////////////////// GEANT4 PROCESS METHODS //////////////
    /////////////////////////////////////////////////////////
public:
    // DO IT
    virtual G4VParticleChange* PostStepDoIt(const G4Track&,
                                            const G4Step& );
    virtual G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                                          const G4Step& aStep){
        return fRegisteredProcess->AtRestDoIt(aTrack,aStep);
    }
    
    // GPIL
    virtual G4double PostStepGetPhysicalInteractionLength (const G4Track&, 
                                                 G4double, G4ForceCondition*);
    virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& aTrack,
                     G4ForceCondition* condition){
        return fRegisteredProcess->AtRestGetPhysicalInteractionLength(aTrack,
               condition);
    };
    
    // GENERAL
    void StartTracking(G4Track*);
    virtual G4bool IsApplicable(const G4ParticleDefinition&);

    // PHYSICS TABLE
    virtual void BuildPhysicsTable(const G4ParticleDefinition&);
    virtual void PreparePhysicsTable(const G4ParticleDefinition&);
    virtual G4bool StorePhysicsTable(const G4ParticleDefinition*,
                                     const G4String&, G4bool);
    virtual G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                                        const G4String&, G4bool);

protected:
    // MFP
    virtual G4double GetMeanFreePath(const G4Track&,
                                     G4double,
                                     G4ForceCondition*);
 
public:
    virtual void EndTracking() {fRegisteredProcess->EndTracking();};
    
    virtual void SetProcessManager(const G4ProcessManager* aPM) 
            {fRegisteredProcess->SetProcessManager(aPM);};
    virtual const G4ProcessManager* GetProcessManager()
            {return fRegisteredProcess->GetProcessManager();};
    virtual void ResetNumberOfInteractionLengthLeft()
            {fRegisteredProcess->ResetNumberOfInteractionLengthLeft();};
    
    virtual void DumpInfo() const {fRegisteredProcess->DumpInfo();};
    
    virtual void SetMasterProcess(G4VProcess* masterP){
        fRegisteredProcess->SetMasterProcess(
           static_cast<XWrapperDiscreteProcess*>(masterP)->fRegisteredProcess
        );
    };
    virtual void BuildWorkerPhysicsTable(const G4ParticleDefinition& aPD)
        {fRegisteredProcess->BuildWorkerPhysicsTable(aPD);};
    virtual void PrepareWorkerPhysicsTable(const G4ParticleDefinition& aPD)
        {fRegisteredProcess->PrepareWorkerPhysicsTable(aPD);};

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

};

#endif

