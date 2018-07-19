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
#include "G4ChannelingOptrChangeCrossSection.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4BOptnChangeCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4VProcess.hh"

#include "Randomize.hh"

#include "G4InteractionLawPhysical.hh"

#include "G4ChannelingTrackData.hh"
#include "G4EmProcessSubType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingOptrChangeCrossSection::G4ChannelingOptrChangeCrossSection(G4String particleName,
                                                                       G4String         name)
:G4VBiasingOperator(name),
fChannelingID(-1),
fSetup(true){
    fParticleToBias = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
    
    if ( fParticleToBias == 0 )
    {
        G4ExceptionDescription ed;
        ed << "Particle `" << particleName << "' not found !" << G4endl;
        G4Exception("G4ChannelingOptrChangeCrossSection(...)",
                    "G4Channeling",
                    JustWarning,
                    ed);
    }
    
    fProcessToDensity["channeling"] = fDensityRatioNone;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingOptrChangeCrossSection::~G4ChannelingOptrChangeCrossSection(){
    for ( std::map< const G4BiasingProcessInterface*, G4BOptnChangeCrossSection* >::iterator
         it = fChangeCrossSectionOperations.begin() ;
         it != fChangeCrossSectionOperations.end() ;
         it++ ) delete (*it).second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChannelingOptrChangeCrossSection::StartRun(){
    if ( fSetup ){
        const G4ProcessManager* processManager = fParticleToBias->GetProcessManager();
        const G4BiasingProcessSharedData* sharedData =
        G4BiasingProcessInterface::GetSharedData( processManager );
        if ( sharedData ){
            for ( size_t i = 0 ; i < (sharedData->GetPhysicsBiasingProcessInterfaces()).size(); i++ ){
                const G4BiasingProcessInterface* wrapperProcess =
                (sharedData->GetPhysicsBiasingProcessInterfaces())[i];
                G4String processName = wrapperProcess->GetWrappedProcess()->GetProcessName();
                G4String operationName = "channelingChangeXS-" + processName;
                fChangeCrossSectionOperations[wrapperProcess] =
                new G4BOptnChangeCrossSection(operationName);
                
                G4ProcessType type = wrapperProcess->GetWrappedProcess()->GetProcessType();
                G4int subType = wrapperProcess->GetWrappedProcess()->GetProcessSubType();

                switch (type) {
                    case fNotDefined:
                        fProcessToDensity[processName] = fDensityRatioNotDefined;
                        break;
                    case fTransportation:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fElectromagnetic:
                        if(subType == fCoulombScattering ||
                           subType == fMultipleScattering){
                            fProcessToDensity[processName] = fDensityRatioNuD;
                        }
                        if(subType == fIonisation ||
                           subType == fPairProdByCharged ||
                           subType == fAnnihilation ||
                           subType == fAnnihilationToMuMu ||
                           subType == fAnnihilationToHadrons){
                            fProcessToDensity[processName] = fDensityRatioElD;
                        }
                        if(subType == fBremsstrahlung ||
                           subType == fNuclearStopping){
                            fProcessToDensity[processName] = fDensityRatioNuDElD;
                        }
                        
                        if(subType == fCerenkov ||
                           subType == fScintillation ||
                           subType == fSynchrotronRadiation ||
                           subType == fTransitionRadiation){
                            fProcessToDensity[processName] = fDensityRatioNone;
                        }
                        if(subType == fRayleigh ||
                           subType == fPhotoElectricEffect ||
                           subType == fComptonScattering ||
                           subType == fGammaConversion ||
                           subType == fGammaConversionToMuMu){
                            fProcessToDensity[processName] = fDensityRatioNone;
                        }
                        break;
                    case fOptical:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fHadronic:
                        fProcessToDensity[processName] = fDensityRatioNuD;
                        break;
                    case fPhotolepton_hadron:
                        fProcessToDensity[processName] = fDensityRatioNuD;
                        break;
                    case fGeneral:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fDecay:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fParameterisation:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fUserDefined:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fParallel:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fPhonon:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    case fUCN:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                    default:
                        fProcessToDensity[processName] = fDensityRatioNone;
                        break;
                }
            }
        }
        fSetup = false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation*
G4ChannelingOptrChangeCrossSection::ProposeOccurenceBiasingOperation(const G4Track*            track,
                                                                     const G4BiasingProcessInterface*
                                                                     callingProcess)
{
    if ( track->GetDefinition() != fParticleToBias ) return 0;
    
    G4double analogInteractionLength =
    callingProcess->GetWrappedProcess()->GetCurrentInteractionLength();
    if ( analogInteractionLength > DBL_MAX/10. ) return 0;
    
    G4double analogXS = 1./analogInteractionLength;
    
    if(fChannelingID==-1){
        fChannelingID = G4PhysicsModelCatalog::GetIndex("channeling");
    }
    G4ChannelingTrackData* trackdata =
    (G4ChannelingTrackData*)(track->GetAuxiliaryTrackInformation(fChannelingID));
    if(trackdata==nullptr) return 0;
    
    G4double XStransformation = 1.;
    auto search = fProcessToDensity.find(callingProcess->GetWrappedProcess()->GetProcessName());
    if(search != fProcessToDensity.end()) {
        switch (search->second) {
            case fDensityRatioNuDElD:
                XStransformation = trackdata->GetDensity();
                break;
            case fDensityRatioNuD:
                XStransformation = trackdata->GetNuD();
                break;
            case fDensityRatioElD:
                XStransformation = trackdata->GetElD();
                break;
            case fDensityRatioNone:
                return 0;
                break;
            case fDensityRatioNotDefined:
                return 0;
                break;
            default:
                return 0;
                break;
        }
    }
    else{
        XStransformation = trackdata->GetDensity();
    }

    G4BOptnChangeCrossSection*   operation = fChangeCrossSectionOperations[callingProcess];
    G4VBiasingOperation* previousOperation = callingProcess->GetPreviousOccurenceBiasingOperation();
    
    if ( previousOperation == 0 ){
        operation->SetBiasedCrossSection( XStransformation * analogXS );
        operation->Sample();
    }
    else{
        if (  previousOperation != operation ){
            G4ExceptionDescription ed;
            ed << " Logic problem in operation handling !" << G4endl;
            G4Exception("G4ChannelingOptrChangeCrossSection::ProposeOccurenceBiasingOperation(...)",
                        "G4Channeling",
                        JustWarning,
                        ed);
            return 0;
        }
        if ( operation->GetInteractionOccured() ){
            operation->SetBiasedCrossSection( XStransformation * analogXS );
            operation->Sample();
        }
        else{
            operation->UpdateForStep( callingProcess->GetPreviousStepSize() );
            operation->SetBiasedCrossSection( XStransformation * analogXS );
            operation->UpdateForStep( 0.0 );
        }
    }
    
    return operation;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChannelingOptrChangeCrossSection::
OperationApplied(const G4BiasingProcessInterface*           callingProcess,
                 G4BiasingAppliedCase,
                 G4VBiasingOperation*             occurenceOperationApplied,
                 G4double,
                 G4VBiasingOperation*,
                 const G4VParticleChange*                                  )
{
    G4BOptnChangeCrossSection* operation = fChangeCrossSectionOperations[callingProcess];
    if ( operation ==  occurenceOperationApplied ) operation->SetInteractionOccured();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
