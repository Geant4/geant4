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
// $Id: G4ITModelProcessor.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITModelProcessor.hh"
#include "G4VITTimeStepper.hh"
#include "G4VITReactionProcess.hh"

std::map<const G4Track*, G4bool> G4ITModelProcessor::fHasReacted ;

G4ITModelProcessor::G4ITModelProcessor()
{
    //ctor
    fpTrack = 0;
    fpModelHandler = 0;
    fpModel = 0;
    fInitialized = false;
    fpModelManager = 0;
    fCurrentModel.assign(G4ITType::size(), std::vector<G4VITModel*>());

    for(int i = 0 ; i < (int) G4ITType::size() ; i++)
    {
        fCurrentModel[i].assign(G4ITType::size(),0);
    }
    fUserMinTimeStep = -1.;
}

G4ITModelProcessor::~G4ITModelProcessor()
{
    //dtor
//    if(fpModelHandler) delete fpModelHandler; deleted by G4ITStepManager
    fCurrentModel.clear();
    fReactionInfo.clear();
}

// Should not be used
G4ITModelProcessor::G4ITModelProcessor(const G4ITModelProcessor& /*other*/)
{
    //copy ctorr
    fpTrack = 0;
    fpModelHandler = 0;
    fpModel = 0;
    fInitialized = false;
    fpModelManager = 0;
    fUserMinTimeStep = -1.;
}

// Should not be used
G4ITModelProcessor& G4ITModelProcessor::operator=(const G4ITModelProcessor& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}
//______________________________________________________________________________
void G4ITModelProcessor::Initialize()
{
    fpModelHandler->Initialize();
    fInitialized = true;
}

//______________________________________________________________________________
void G4ITModelProcessor::InitializeStepper(const G4double& currentGlobalTime,
                                           const G4double& userMinTime)
{
    // G4cout << "G4ITModelProcessor::InitializeStepper" << G4endl;
    if(fpModelHandler==0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "No G4ITModelHandler was passed to the modelProcessor.";
        G4Exception("G4ITModelProcessor::InitializeStepper","ITModelProcessor002",
                    FatalErrorInArgument,exceptionDescription);
    }
    const std::vector<std::vector<G4ITModelManager*> >* modelManager = fpModelHandler
            ->GetAllModelManager();

    if(modelManager==0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "No G4ITModelManager was register to G4ITModelHandler.";
        G4Exception("G4ITModelProcessor::InitializeStepper","ITModelProcessor003",
                    FatalErrorInArgument,exceptionDescription);
    }

    int nbModels1 = modelManager->size() ;

    G4VITTimeStepper::SetTimes(currentGlobalTime, userMinTime) ;

    // TODO !!!
    //    if( nbModels1 != 1 || (nbModels1 == 1 && !fpModelManager) )
    {
        int nbModels2 = -1;
        G4VITModel* model = 0;
        G4ITModelManager* modman = 0;

        for(int i = 0 ; i < nbModels1 ; i++)
        {
            nbModels2 = (*modelManager)[i].size();

            for(int j = 0 ; j <= i ; j++)
            {
                modman = (*modelManager)[i][j];

                if(modman == 0) continue ;

                model       =  modman -> GetModel(currentGlobalTime);
                G4VITTimeStepper* stepper   = model->GetTimeStepper() ;

//                stepper -> PrepareForAllProcessors() ;
                stepper -> Prepare() ;
                fCurrentModel[i][j] = model;
            }
        }

        if(nbModels1 == 1 && nbModels2 ==1)
        {
            fpModelManager = modman;
            fpModel = model;
        }
        else fpModel = 0;
    }
}

//______________________________________________________________________________
void G4ITModelProcessor::CalculateTimeStep(const G4Track* track, const G4double userMinTimeStep)
{
    // G4cout  << "G4ITModelProcessor::CalculateStep" << G4endl;
    CleanProcessor();
    if(track == 0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "No track was passed to the method (track == 0).";
        G4Exception("G4ITModelProcessor::CalculateStep","ITModelProcessor004",
                    FatalErrorInArgument,exceptionDescription);
    }
    SetTrack(track);
    fUserMinTimeStep = userMinTimeStep ;

    DoCalculateStep();
}

//______________________________________________________________________________
void G4ITModelProcessor::DoCalculateStep()
{
    if(fpModel) // ie only one model has been declared and will be used
    {
        fpModel -> GetTimeStepper()->CalculateStep(*fpTrack, fUserMinTimeStep);
    }
    else // ie many models have been declared and will be used
    {
        std::vector<G4VITModel*>& model = fCurrentModel[GetIT(fpTrack)->GetITType()];

        for(int i =0 ; i < (int) model.size() ; i++)
        {
            if(model[i] == 0) continue;
            model[i]->GetTimeStepper()->CalculateStep(*fpTrack, fUserMinTimeStep);
        }
    }
}

//______________________________________________________________________________
void G4ITModelProcessor::FindReaction(std::map<G4Track*, G4TrackVectorHandle>* tracks,
                                      const double currentStepTime,
                                      const double previousStepTime,
                                      const bool reachedUserStepTimeLimit)
{
    // DEBUG
    //    G4cout << "G4ITReactionManager::FindReaction" << G4endl;
    if(tracks == 0)       return ;

    if(fpModelHandler->GetAllModelManager()->empty()) return ;

    std::map<G4Track*, G4TrackVectorHandle>::iterator tracks_i = tracks->begin();;

    for(tracks_i = tracks->begin() ; tracks_i != tracks-> end() ; tracks_i ++)
    {
        /// Get track A
        G4Track* trackA = tracks_i->first;

        if(trackA == 0)         continue;

        std::map<const G4Track*, G4bool>::iterator it_hasReacted = fHasReacted.find(trackA);
        if(it_hasReacted != fHasReacted.end()) continue;
        if(trackA->GetTrackStatus() == fStopAndKill) continue;

        G4IT* ITA = GetIT(trackA);
        G4ITType ITypeA = ITA -> GetITType();

        const std::vector<G4VITModel*> model = fCurrentModel[ITypeA];

        G4TrackVectorHandle& trackB_vector = tracks_i->second ;
        std::vector<G4Track*>::iterator trackB_i = trackB_vector->begin();

        G4Track* trackB = 0 ;
        G4ITType ITypeB(-1);
        G4VITReactionProcess* process = 0;
        G4ITReactionChange* changes = 0;

        for(; trackB_i != trackB_vector->end() ; trackB_i++)
        {
            trackB = *trackB_i;

            if(trackB == 0)         continue;
            it_hasReacted = fHasReacted.find(trackB);
            if(it_hasReacted != fHasReacted.end()) continue;
            if(trackB->GetTrackStatus() == fStopAndKill) continue;

            // DEBUG
            //             G4cout << "Couple : " << trackA->GetParticleDefinition->GetParticleName() << " ("
            //                        << trackA->GetTrackID() << ")   "
            //                        << trackB->GetParticleDefinition->GetParticleName() << " ("
            //                        << trackB->GetTrackID() << ")"
            //                        << G4endl;

            if(trackB == trackA)
            {
                G4ExceptionDescription exceptionDescription ;
                exceptionDescription << "The IT reaction process sent back a reaction between trackA and trackB. ";
                exceptionDescription << "The problem is trackA == trackB";
                G4Exception("G4ITModelProcessor::FindReaction","ITModelProcessor005",
                            FatalErrorInArgument,exceptionDescription);
            }

            G4IT* ITB = GetIT(trackB);
            G4ITType ITypeBtmp = ITB -> GetITType();

            if(ITypeB != ITypeBtmp)
            {
                ITypeB = ITypeBtmp ;

                if(model[ITypeB])
                    process = model[ITypeB]->GetReactionProcess();
            }

            if(process && process -> TestReactibility(*trackA, *trackB,
                                                      currentStepTime, previousStepTime,
                                                      reachedUserStepTimeLimit))
            {
                changes = process->MakeReaction(*trackA, *trackB);
            }

            if(changes)
            {
                fHasReacted[trackA] = true;
                fHasReacted[trackB] = true;
                changes -> GetTrackA();
                changes -> GetTrackB();

                fReactionInfo.push_back(changes);

                process->ResetChanges();
                changes = 0;

                break;
            }
        }
    }

    fHasReacted.clear();
}
