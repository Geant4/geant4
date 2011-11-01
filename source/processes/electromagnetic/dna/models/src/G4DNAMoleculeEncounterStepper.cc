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
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fReactionTable)),
    fVerbose(0)
{}

G4DNAMoleculeEncounterStepper::~G4DNAMoleculeEncounterStepper()
{}

G4DNAMoleculeEncounterStepper::G4DNAMoleculeEncounterStepper(const G4DNAMoleculeEncounterStepper& right) :
    G4VITTimeStepper(right),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fReactionTable))
{
    fVerbose                 = right.fVerbose ;
    fMolecularReactionTable  = right.fMolecularReactionTable;
}

void G4DNAMoleculeEncounterStepper::PrepareForAllProcessors()
{
//    G4cout << "G4DNAMoleculeEncounterStepper::PrepareForAllProcessors" << G4endl;
    G4ITManager<G4Molecule>::Instance()->UpdatePositionMap();
}

G4double G4DNAMoleculeEncounterStepper::CalculateStep(const G4Track& trackA, const G4double& userMinTimeStep)
{
//    G4cout << "G4MoleculeEncounterStepper::CalculateStep, time :" << G4ITStepManager::Instance()->GetGlobalTime()  << G4endl;
    fUserMinTimeStep = userMinTimeStep ;
    if(fReactants) fReactants = 0 ;

    fReactants = new vector<G4Track*>();

    G4Molecule* moleculeA = GetMolecule(trackA);

    //__________________________________________________________________
    // Retrieve general informations for making reactions
    const vector<const G4Molecule*>* reactivesVector =
            fMolecularReactionTable -> CanReactWith(moleculeA);

    if(!reactivesVector)
    {
        // G4cout << "G4MoleculeEncounterStepper::CalculateStep will return DBL_MAX" << G4endl;
        return DBL_MAX;
    }

    fReactionModel -> Initialise(moleculeA, trackA) ;
    G4int nbReactives = reactivesVector->size();

    if(nbReactives == 0)
    {
        G4cout << "!!!!!!!!!!!!!!!!!!!!"<<G4endl;
        G4cout << "!!! WARNING" << G4endl;
        G4cout << "Return for molecule " << moleculeA -> GetName ()  << G4endl ;
        G4cout << "G4MoleculeEncounterStepper::CalculateStep will return DBL_MAX" << G4endl;
        G4cout << "because "<< moleculeA -> GetName () << " does not have any reactants "
               << "given in the reaction table." << G4endl;
        G4cout << "!!!!!!!!!!!!!!!!!!!!"<<G4endl;
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

    fSampledMinTimeStep = DBL_MAX;
    G4bool hasAlreadyReachNullTime = false;

    if(fVerbose>=4)
    {
        G4cout << "___________________________________" << G4endl;
        G4cout << "Incident molecule : " << moleculeA->GetName() << G4endl;
    }
    //__________________________________________________________________
    // Start looping on possible reactants
    for (G4int i=0 ; i<nbReactives ; i++)
    {
        const G4Molecule* moleculeB = (*reactivesVector)[i];

        G4double DA = moleculeA->GetDiffusionCoefficient() ;
        G4double DB = moleculeB->GetDiffusionCoefficient() ;


        //______________________________________________________________
        // Retrieve reaction range
        G4double R = -1 ; // reaction Range
        R = fReactionModel -> GetReactionRadius(i);

        //______________________________________________________________
        // Use KdTree algorithm to find closest reactants
        G4KDTreeResultHandle results (G4ITManager<G4Molecule>::Instance()
                                      -> FindNearest(moleculeA, moleculeB));

        if(!results)
        {
            // DEBUG
            //                G4cout << " Pas de molecule trouvée : " << moleculeB->GetName() << G4endl;
            continue ;
        }

        for(results->Rewind();
            !results->End();
            results->Next())
        {

            G4IT* reactiveB = (G4IT*) results->GetItemData() ;
            if (reactiveB->GetTrack()->GetTrackStatus() != fAlive)
            {
                G4ExceptionDescription exceptionDescription ("The track status of one of the nearby reactants is not fAlive");
                G4Exception("G4DNAMoleculeEncounterStepper::CalculateStep","MoleculeEncounterStepper001",
                            FatalErrorInArgument,exceptionDescription);
                continue ;
            }

            if (reactiveB==0)
            {
                //  DEBUG
                //  G4cout<<"Continue 1"<<G4endl;
                continue ;
            }

            G4Track *trackB = reactiveB->GetTrack();

            if(trackB == &trackA)
            {
                // DEBUG
                G4cerr << "Molecule B is of type : "
                       << moleculeB->GetName()
                       << " trackID : "
                       << trackB->GetTrackID()<< G4endl;
                G4cerr << "Molecule A is of type : "
                       << moleculeA->GetName()
                       << " trackID : "
                       << trackA.GetTrackID() << G4endl;

                G4ExceptionDescription exceptionDescription ("The track you are requested nearby reactants for ");
                exceptionDescription << " and the nearby reactant returned by the ITManager are the same" ;
                G4Exception("G4DNAMoleculeEncounterStepper::CalculateStep","MoleculeEncounterStepper002",
                            FatalErrorInArgument,exceptionDescription);

            }

            if(fabs(trackB->GetGlobalTime() - trackA.GetGlobalTime()) > trackA.GetGlobalTime()*(1-1/100) )
            {
                // DEBUG
                G4cerr<< "trackB->GetGlobalTime() != trackA.GetGlobalTime()"
                      << G4endl;

                G4cerr << "trackID : " << trackA.GetTrackID()
                       << "\t Name :" << moleculeA->GetName()
                       <<"\t trackA->GetGlobalTime() = "
                       << G4BestUnit(trackA.GetGlobalTime(), "Time") << G4endl;

                G4cerr << "trackID : " << trackB->GetTrackID()
                       << "\t Name :" << moleculeB->GetName()
                       << "\t trackB->GetGlobalTime() = "
                       << G4BestUnit(trackB->GetGlobalTime(), "Time")<< G4endl;

                G4ExceptionDescription exceptionDescription ("The tracks are not synchronized in time");
                G4Exception("G4DNAMoleculeEncounterStepper::CalculateStep","MoleculeEncounterStepper003",
                            FatalErrorInArgument,exceptionDescription);
            }

            G4double r2 = results->GetDistanceSqr() ;

            if(fVerbose>=4)
            {
                G4cout <<"Reaction : Interaction Range = "
                      << G4BestUnit(R, "Length")<<G4endl;
                G4cout <<"Real distance between reactants  = "
                       << G4BestUnit((trackA.GetPosition() - trackB->GetPosition()).mag(), "Length")<<G4endl;
                G4cout <<"Distance between reactants calculated by nearest neighbor algorithm = "
                       << G4BestUnit(sqrt(r2), "Length")<<G4endl;
            }

            if(r2 <= R*R)
            {
                if(hasAlreadyReachNullTime == false)
                {
                    fReactants->clear();
                    hasAlreadyReachNullTime = true;
                }
                fReactants->push_back(trackB);
                fSampledMinTimeStep = 0.;

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

    //    DEBUG
//    G4cout << "G4MoleculeEncounterStepper::CalculateStep will return :"
//           << G4BestUnit(fSampledMinTimeStep, "Time")<< G4endl;
    return fSampledMinTimeStep ;
}
