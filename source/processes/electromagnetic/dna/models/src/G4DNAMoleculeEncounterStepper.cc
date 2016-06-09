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
// $Id: G4DNAMoleculeEncounterStepper.cc 64374 2012-10-31 16:37:23Z gcosmo $
//
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4VDNAReactionModel.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4H2O.hh"
#include "G4UnitsTable.hh"

using namespace std;

G4DNAMoleculeEncounterStepper::G4DNAMoleculeEncounterStepper() :
    G4VITTimeStepper(),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable)),
    fReactionModel(0)
{
    fVerbose = 0;
    fHasAlreadyReachedNullTime = false;
}

G4DNAMoleculeEncounterStepper& G4DNAMoleculeEncounterStepper::operator=(const G4DNAMoleculeEncounterStepper& rhs)
{
    if(this == &rhs) return *this;
    fReactionModel              = 0;
    fVerbose                    = rhs.fVerbose;
    fMolecularReactionTable     = rhs.fMolecularReactionTable;
    fHasAlreadyReachedNullTime  = false;
    return *this;
}

G4DNAMoleculeEncounterStepper::~G4DNAMoleculeEncounterStepper()
{}

G4DNAMoleculeEncounterStepper::G4DNAMoleculeEncounterStepper(const G4DNAMoleculeEncounterStepper& right) :
    G4VITTimeStepper(right),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
    fVerbose                    = right.fVerbose ;
    fMolecularReactionTable     = right.fMolecularReactionTable;
    fReactionModel              = 0;
    fHasAlreadyReachedNullTime  = false;
}

void G4DNAMoleculeEncounterStepper::Prepare()
{
    // DEBUG
    //    G4cout << "G4DNAMoleculeEncounterStepper::PrepareForAllProcessors" << G4endl;
    G4VITTimeStepper::Prepare();
    G4ITManager<G4Molecule>::Instance()->UpdatePositionMap();
}

G4double G4DNAMoleculeEncounterStepper::CalculateStep(const G4Track& trackA, const G4double& userMinTimeStep)
{
    // DEBUG
    //    G4cout << "G4MoleculeEncounterStepper::CalculateStep, time :" << G4ITTrackHolder::Instance()->GetGlobalTime()  << G4endl;

    G4Molecule* moleculeA = GetMolecule(trackA);


#ifdef G4VERBOSE
    if(fVerbose)
    {
        G4cout << "_______________________________________________________________________" << G4endl;
        G4cout << "G4DNAMoleculeEncounterStepper::CalculateStep" << G4endl;
        G4cout << "Incident molecule : " << moleculeA->GetName()
               << " (" << trackA.GetTrackID() << ") "
               << G4endl;
    }
#endif

    //__________________________________________________________________
    // Retrieve general informations for making reactions
    const vector<const G4Molecule*>* reactivesVector =
            fMolecularReactionTable -> CanReactWith(moleculeA);

    if(!reactivesVector)
    {
#ifdef G4VERBOSE
        //    DEBUG
        if(fVerbose > 1)
        {
            G4cout << "!!!!!!!!!!!!!!!!!!!!"<<G4endl;
            G4cout << "!!! WARNING" << G4endl;
            G4cout << "G4MoleculeEncounterStepper::CalculateStep will return infinity for the reaction because the molecule "
                   << moleculeA->GetName()
                   << " does not have any reactants given in the reaction table."
                   << G4endl;
            G4cout << "!!!!!!!!!!!!!!!!!!!!"<<G4endl;
        }
#endif
        return DBL_MAX;
    }

    G4int nbReactives = reactivesVector->size();

    if(nbReactives == 0)
    {
#ifdef G4VERBOSE
        //    DEBUG
        if(fVerbose)
        {
            G4cout << "!!!!!!!!!!!!!!!!!!!!"<<G4endl;
            G4cout << "!!! WARNING" << G4endl;
            G4cout << "G4MoleculeEncounterStepper::CalculateStep will return infinity for the reaction because the molecule "
                   << moleculeA->GetName()
                   << " does not have any reactants given in the reaction table."
                   << "This message can also result from a wrong implementation of the reaction table."
                   << G4endl;
            G4cout << "!!!!!!!!!!!!!!!!!!!!"<<G4endl;
        }
#endif
        return DBL_MAX;
    }
    // DEBUG
    //    else
    //    {
    //        G4cout << "nb reactants : " << nbReactives << " pour mol "<< moleculeA -> GetName () << G4endl;
    //        for(int k=0 ; k < nbReactives ; k++)
    //        {
    //            G4cout << (*reactivesVector)[k]->GetName() << G4endl;
    //        }
    //    }

    fUserMinTimeStep = userMinTimeStep ;
    if(fReactants) fReactants = 0 ;
    fReactants = new vector<G4Track*>();

    fSampledMinTimeStep = DBL_MAX;
    fHasAlreadyReachedNullTime = false;

    fReactionModel -> Initialise(moleculeA, trackA) ;

    //__________________________________________________________________
    // Start looping on possible reactants
    for (G4int i=0 ; i<nbReactives ; i++)
    {
        const G4Molecule* moleculeB = (*reactivesVector)[i];

        //______________________________________________________________
        // Retrieve reaction range
        G4double R = -1 ; // reaction Range
        R = fReactionModel -> GetReactionRadius(i);

        //______________________________________________________________
        // Use KdTree algorithm to find closest reactants
        G4KDTreeResultHandle results (G4ITManager<G4Molecule>::Instance()
                                      -> FindNearest(moleculeA, moleculeB));

        RetrieveResults(trackA,moleculeA,moleculeB,R,results, true);
    }

#ifdef G4VERBOSE
    //    DEBUG
    if(fVerbose)
    {
        G4cout << "G4MoleculeEncounterStepper::CalculateStep will finally return :"
               << G4BestUnit(fSampledMinTimeStep, "Time")<< G4endl;

        if(fVerbose > 1)
        {
            G4cout << "TrackA: " << moleculeA->GetName() << " (" << trackA.GetTrackID() << ") can react with: ";

            vector<G4Track*>::iterator it;
            for(it = fReactants->begin() ; it != fReactants->end() ; it++)
            {
                G4Track* trackB = *it;
                G4cout << GetMolecule(trackB)->GetName() << " ("
                       << trackB->GetTrackID() << ") \t ";
            }
            G4cout << G4endl;
        }
    }
#endif
    return fSampledMinTimeStep ;
}


void G4DNAMoleculeEncounterStepper::RetrieveResults(const G4Track& trackA, const G4Molecule* moleculeA,
                                                    const G4Molecule* moleculeB,  const G4double R,
                                                    G4KDTreeResultHandle& results, G4bool iterate)
{
    if(!results)
    {
#ifdef G4VERBOSE
        // DEBUG
        if(fVerbose > 1)
        {
            G4cout << "No molecule " << moleculeB->GetName()
                   << " found to react with "
                   << moleculeA->GetName()
                   << G4endl;
        }
#endif
        return ;
    }

    G4double DA = moleculeA->GetDiffusionCoefficient() ;
    G4double DB = moleculeB->GetDiffusionCoefficient() ;

    for(results->Rewind();
        !results->End();
        results->Next())
    {

        G4IT* reactiveB = (G4IT*) results->GetItemData() ;

        if (reactiveB==0)
        {
            //  DEBUG
            //  G4cout<<"Continue 1"<<G4endl;
            continue ;
        }

        G4Track *trackB = reactiveB->GetTrack();

        if(trackB == 0)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "The reactant B found using the ITManager does not have a valid track "
                                 << " attached to it. If this is done on purpose, please do not record this "
                                 << " molecule in the ITManager."
                                 << G4endl;
            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults","MoleculeEncounterStepper001",
                        FatalErrorInArgument,exceptionDescription);
            continue ;
        }

        if (trackB->GetTrackStatus() != fAlive)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "The track status of one of the nearby reactants is not fAlive" << G4endl;
            exceptionDescription << "The incomming trackID "
                                 << "(trackA entering in G4DNAMoleculeEncounterStepper and "
                                 << "for which you are looking reactant for) is : "
                                 << trackA.GetTrackID()  <<"("<< GetMolecule(trackA)->GetName()<<")" << G4endl;
            exceptionDescription << "And the trackID of the reactant (trackB) is: "
                                 << trackB->GetTrackID()  <<"("<< GetMolecule(trackB)->GetName()<<")" << G4endl;
            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults","MoleculeEncounterStepper002",
                        FatalErrorInArgument,exceptionDescription);
            continue ;
        }

        if(trackB == &trackA)
        {
            // DEBUG
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "A track is reacting with itself (which is impossible) ie trackA == trackB"
                                 << G4endl ;
            exceptionDescription << "Molecule A (and B) is of type : "
                                 << moleculeA->GetName()
                                 << " with trackID : "
                                 << trackA.GetTrackID() << G4endl;

            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults","MoleculeEncounterStepper003",
                        FatalErrorInArgument,exceptionDescription);

        }

        if(fabs(trackB->GetGlobalTime() - trackA.GetGlobalTime()) > trackA.GetGlobalTime()*(1-1/100) )
        {
            // DEBUG
            G4ExceptionDescription exceptionDescription;
            exceptionDescription << "The interacting tracks are not synchronized in time"<< G4endl;
            exceptionDescription<< "trackB->GetGlobalTime() != trackA.GetGlobalTime()"
                                << G4endl;

            exceptionDescription << "trackA : trackID : " << trackA.GetTrackID()
                                 << "\t Name :" << moleculeA->GetName()
                                 <<"\t trackA->GetGlobalTime() = "
                                << G4BestUnit(trackA.GetGlobalTime(), "Time") << G4endl;

            exceptionDescription << "trackB : trackID : " << trackB->GetTrackID()
                                 << "\t Name :" << moleculeB->GetName()
                                 << "\t trackB->GetGlobalTime() = "
                                 << G4BestUnit(trackB->GetGlobalTime(), "Time")<< G4endl;

            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults","MoleculeEncounterStepper004",
                        FatalErrorInArgument,exceptionDescription);
        }

        G4double r2 = results->GetDistanceSqr() ;
#ifdef G4VERBOSE
        if(fVerbose > 1)
        {
            G4cout << "\t ************************************************** " << G4endl;
            G4cout <<"\t Reaction between "
                  << moleculeA->GetName() << " (" << trackA.GetTrackID() << ") "
                   << " & " << moleculeB->GetName() << " (" << trackB->GetTrackID() << "), "
                   << "Interaction Range = "
                   << G4BestUnit(R, "Length")<<G4endl;
            G4cout <<"\t Real distance between reactants  = "
                  << G4BestUnit((trackA.GetPosition() - trackB->GetPosition()).mag(), "Length")<<G4endl;
            G4cout <<"\t Distance between reactants calculated by nearest neighbor algorithm = "
                  << G4BestUnit(sqrt(r2), "Length")<<G4endl;
//            G4cout << " ***** " << G4endl;
        }
#endif

        if(r2 <= R*R)
        {
            // Entering in this condition may due to the fact that molecules are very close
            // to each other
            // Therefore, if we consider only the nearby reactant, this one may have already
            // react. Instead, we will take all possible reactants that satisfy the condition r<R

            if(fHasAlreadyReachedNullTime == false)
            {
                fReactants->clear();
                fHasAlreadyReachedNullTime = true;
            }

            if(iterate) // First call (call from outside this method)
            {
                fSampledMinTimeStep = 0.;
                G4KDTreeResultHandle results2 (G4ITManager<G4Molecule>::Instance()
                                          -> FindNearestInRange(moleculeA, moleculeB,R));
                RetrieveResults(trackA, moleculeA, moleculeB, R, results2, false);
            }
            else // Other calls (call from this method)
            {
                fReactants->push_back(trackB);
            }

            continue;
        }
        else
        {
            G4double r = sqrt(r2);
            G4double tempMinET = pow(r - R,2)
                    /(16 * (DA + DB + 2*sqrt(DA*DB)));

            if(tempMinET <= fSampledMinTimeStep)
            {
                if(tempMinET <= fUserMinTimeStep)
                {
                    if(fSampledMinTimeStep > fUserMinTimeStep)
                        fReactants->clear();
                    fSampledMinTimeStep = fUserMinTimeStep;
                    fReactants->push_back(trackB);

                }
                else
                {
                    fSampledMinTimeStep = tempMinET;
                    if(tempMinET < fSampledMinTimeStep)
                        fReactants->clear();
                    fReactants->push_back(trackB);
                }
            }
        }
    }
}
