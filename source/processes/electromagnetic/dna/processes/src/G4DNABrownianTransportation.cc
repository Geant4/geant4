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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNABrownianTransportation.hh"
#include "G4Molecule.hh"
#include "G4RandomDirection.hh"
#include "G4ParticleTable.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

using namespace std;

#ifndef State
#define State(theXInfo) (fTransportationState->theXInfo)
#endif

//COLOR FOR DEBUGING
//#define RED_ON_WHITE  "\033[0;31m"
//#define GREEN "\033[32;40m"
#define GREEN_ON_BLUE "\033[1;32;44m"
#define RESET "\033[0m"

G4DNABrownianTransportation::G4DNABrownianTransportation(const G4String& aName, G4int verbosity) :
    G4ITTransportation(aName, verbosity)
{
    //ctor
    fpSafetyHelper = G4TransportationManager::GetTransportationManager()
            ->GetSafetyHelper();
    SetProcessSubType(61);
    verboseLevel = 1;
    fNistWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
}

G4DNABrownianTransportation::~G4DNABrownianTransportation()
{;}

G4DNABrownianTransportation::G4DNABrownianTransportation(const G4DNABrownianTransportation& right) :
    G4ITTransportation(right)
{
    //copy ctor
    fpSafetyHelper = G4TransportationManager::GetTransportationManager()
            ->GetSafetyHelper();
}

G4DNABrownianTransportation& G4DNABrownianTransportation::operator=(const G4DNABrownianTransportation& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void G4DNABrownianTransportation::BuildPhysicsTable(const G4ParticleDefinition& particle)
{
    fpSafetyHelper->InitialiseHelper();

    if(verboseLevel > 0)
    {
        G4cout << G4endl << GetProcessName() << ":   for  "
               << setw(24) << particle.GetParticleName()
               << "\tSubType= " << GetProcessSubType()   << G4endl;
    }
}

void G4DNABrownianTransportation::ComputeStep(const G4Track& track,
                                              const G4Step& /*step*/,
                                              const double timeStep,
                                              double& spaceStep)
{
    // G4cout << "G4ITBrownianTransportation::ComputeStep" << G4endl;

    // If this method is called, this step
    // cannot be geometry limited.
    // In case the step is limited by the geometry,
    // this method should not be called.
    // The fTransportEndPosition calculated in
    // the method AlongStepIL should be taken
    // into account.
    // In order to do so, the flag IsLeadingStep
    // is on. Meaning : this track has the minimum
    // interaction length over all others.
    if(GetIT(track)->GetTrackingInfo()->IsLeadingStep())
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "ComputeStep is called while the track has the minimum interaction time";
        exceptionDescription << " so it should not recompute a timeStep ";
        G4Exception("G4DNABrownianTransportation::ComputeStep","G4DNABrownianTransportation001",
                    FatalErrorInArgument,exceptionDescription);
    }

    State(fGeometryLimitedStep) = false;
    // TODO : generalize this process to all kind of brownian objects
    //    G4ITBrownianObject* ITBrown = GetITBrownianObject(track) ;
    //    G4double diffCoeff = ITBrown->GetDiffusionCoefficient(track.GetMaterial());
    G4Molecule* molecule = GetMolecule(track) ;
    G4double diffCoeff = molecule->GetDiffusionCoefficient();

    if(timeStep > 0)
    {
        G4double x = G4RandGauss::shoot(0,sqrt(2*diffCoeff*timeStep));
        G4double y = G4RandGauss::shoot(0,sqrt(2*diffCoeff*timeStep));
        G4double z = G4RandGauss::shoot(0,sqrt(2*diffCoeff*timeStep));

        spaceStep = sqrt(x*x + y*y + z*z);

        State(fTransportEndPosition).set(x + track.GetPosition().x(),
                                         y + track.GetPosition().y(),
                                         z + track.GetPosition().z());
        State(fCandidateEndGlobalTime) = track.GetStep()->GetPreStepPoint()->GetGlobalTime() + timeStep ;
    }
    else
    {
        spaceStep = 0. ;
        State(fTransportEndPosition) =  track.GetPosition() ;
        State(fCandidateEndGlobalTime) = track.GetStep()->GetPreStepPoint()->GetGlobalTime() ;
    }

    State(fEndGlobalTimeComputed) = true ;

#ifdef G4VERBOSE
    //    DEBUG
    if(fVerboseLevel>1)
    {
        G4cout<< GREEN_ON_BLUE
              << "G4ITBrownianTransportation::ComputeStep() : "
              << " trackID : "         << track.GetTrackID()
              << " : Molecule name: "  << molecule-> GetName()
              << G4endl
              << "Diffusion length : " << G4BestUnit(spaceStep, "Length")
              << " within time step : "  << G4BestUnit(timeStep,"Time")
              << RESET
              << G4endl<< G4endl;
    }
#endif
}

G4VParticleChange* G4DNABrownianTransportation::PostStepDoIt( const G4Track& track, const G4Step& step)
{
    G4ITTransportation::PostStepDoIt(track,step);
    Diffusion(track);

#ifdef G4VERBOSE
    //    DEBUG
    if(fVerboseLevel>1)
    {
        G4cout<< GREEN_ON_BLUE
              << "G4ITBrownianTransportation::PostStepDoIt() :"
              << " trackID : "        << track.GetTrackID()
              << " Molecule name: "   << GetMolecule(track)-> GetName() << G4endl;
        G4cout<< "Diffusion length : "<< G4BestUnit(step.GetStepLength(),"Length") <<" within time step : "
              << G4BestUnit(step.GetDeltaTime(),"Time") << "\t"
              << " Current global time : " << G4BestUnit(track.GetGlobalTime(),"Time")
              << RESET
              << G4endl<< G4endl;
    }
#endif
    return &fParticleChange ;
}

void G4DNABrownianTransportation::Diffusion(
    const G4Track& track)
{

#ifdef G4VERBOSE
    // DEBUG
    if (fVerboseLevel>1)
    {
        G4cout<< GREEN_ON_BLUE
              << setw(18)<< "G4DNABrownianTransportation::Diffusion :"
              << setw(8) <<  GetIT(track)->GetName()
              << "\t trackID:"      << track.GetTrackID() <<"\t"
              << " Global Time = "  << G4BestUnit(track.GetGlobalTime(),"Time")
              << RESET
              << G4endl<< G4endl;
    }
#endif

    G4Material* material = track.GetMaterial();
    if (material != fNistWater && material->GetBaseMaterial() != fNistWater)
    {


        G4cout << (track.GetLocalTime()) /s<<G4endl;
        G4cout<< "Step Number :" << track.GetCurrentStepNumber() <<G4endl;

        fParticleChange.ProposeEnergy(0.) ;
        fParticleChange.ProposeTrackStatus(fStopButAlive);

        // Got pb with :
        // fParticleChange.ProposeTrackStatus(fStopAndKill);
        // It makes the tracks going straight without killing them

        return ; // &fParticleChange is the final returned object
    }

    G4double costheta = (2*G4UniformRand()-1);
    G4double theta = acos (costheta);
    G4double phi = 2*pi*G4UniformRand();

    G4double xMomentum = cos(phi)* sin(theta);
    G4double yMomentum = sin(theta)*sin(phi);
    G4double zMomentum = costheta;

    fParticleChange.ProposeMomentumDirection(xMomentum, yMomentum, zMomentum);

    // Alternative
    //fParticleChange.ProposeMomentumDirection(G4RandomDirection());

    return; // &fParticleChange is the final returned object
}


G4double G4DNABrownianTransportation::AlongStepGetPhysicalInteractionLength(
    const G4Track& track,
    G4double , //   previousStepSize
    G4double currentMinimumStep,
    G4double& currentSafety,
    G4GPILSelection* selection)
{
    G4double geometryStepLength(-1), newSafety(-1) ;

    // Initial actions moved to  StartTrack()
    // --------------------------------------
    // Note: in case another process changes touchable handle
    //    it will be necessary to add here (for all steps)
    State(fCurrentTouchableHandle) = track.GetTouchableHandle();

    // GPILSelection is set to default value of CandidateForSelection
    // It is a return value
    //
    *selection = CandidateForSelection ;

    // Get initial Energy/Momentum of the track
    //
    const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
    G4ThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
    G4ThreeVector startPosition          = track.GetPosition() ;

    // The Step Point safety can be limited by other geometries and/or the
    // assumptions of any process - it's not always the geometrical safety.
    // We calculate the starting point's isotropic safety here.
    //
    G4ThreeVector OriginShift = startPosition - State(fPreviousSftOrigin) ;
    G4double      MagSqShift  = OriginShift.mag2() ;
    if( MagSqShift >= sqr(State(fPreviousSafety)) )
    {
        currentSafety = 0.0 ;
    }
    else
    {
        currentSafety = State(fPreviousSafety) - std::sqrt(MagSqShift) ;
    }

    //G4double linearStepLength ;
    if( fShortStepOptimisation && (currentMinimumStep <= currentSafety) )
    {
        // The Step is guaranteed to be taken
        //
        geometryStepLength   = currentMinimumStep ;
        State(fGeometryLimitedStep) = false ;
    }
    else
    {
        newSafety = fpSafetyHelper->ComputeSafety(startPosition);

        // Remember last safety origin & value.
        //
        State(fPreviousSftOrigin) = startPosition ;
        State(fPreviousSafety)    = newSafety ;
        // fpSafetyHelper->SetCurrentSafety( newSafety, startPosition);

        // The safety at the initial point has been re-calculated:
        //
        currentSafety = newSafety ;
    }

    if(currentSafety <= currentMinimumStep)
    {
        geometryStepLength = currentSafety;

        // Remember : brownian motion # straight trajectories
        State(fGeometryLimitedStep) = true ;
        //        State(endpointDistance) = currentSafety ;

        // Calculate final position
        //
        State(fTransportEndPosition) = startPosition+currentSafety*startMomentumDir ;

        // Find out the travelling time from the calculated distance
        // => méthode inverse : d_99
        G4double diffusionCoefficient = GetMolecule(track)
                ->GetDiffusionCoefficient();
        //    G4double diffusionCoefficient = GetITBrownianObject(track)
        //                                    ->GetDiffusionCoefficient(track.GetMaterial());

        // 99 % of the space step distribution is lower than
        // d_99 = 8 * sqrt(D*t)
        // where t is the corresponding time step
        // so by inversion :
        State(theInteractionTimeLeft)     = (currentSafety*currentSafety)/(64 * diffusionCoefficient);
    }
    /*
    else
    {
        //        // TODO !!!
        //        State(endpointDistance) = geometryStepLength ;

        //        // Calculate final position
        //        //
        //        State(fTransportEndPosition) = startPosition+geometryStepLength*startMomentumDir ;
        //        // since the brownian motion is not straight the above formula is not correct
        //        // However, at present time, in DNA no other processes are used in combinaison
        //        // with the brownian motion, so this is not a problem, because the calculated
        //        // position will not be used
    }*/

    // Momentum direction, energy and polarisation are unchanged by transport
    //
    State(fTransportEndMomentumDir)   = startMomentumDir ;
    State(fTransportEndKineticEnergy) = track.GetKineticEnergy() ;
    State(fTransportEndSpin)          = track.GetPolarization();
    State(fParticleIsLooping)         = false ;
    State(fMomentumChanged)           = false ;
    State(fEndGlobalTimeComputed)     = false ;

    // If we are asked to go a step length of 0, and we are on a boundary
    // then a boundary will also limit the step -> we must flag this.
    //
    if( currentMinimumStep == 0.0 )
    {
        if( currentSafety == 0.0 )  State(fGeometryLimitedStep) = true ;
    }

    return geometryStepLength ;
}

//////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)
G4VParticleChange* G4DNABrownianTransportation::AlongStepDoIt( const G4Track& track,
                                                               const G4Step&  step )
{
    static G4int noCalls=0;
    noCalls++;

    fParticleChange.Initialize(track) ;

    fParticleChange.ProposeMomentumDirection(State(fTransportEndMomentumDir)) ;
    fParticleChange.ProposeEnergy(State(fTransportEndKineticEnergy)) ;
    fParticleChange.SetMomentumChanged(State(fMomentumChanged)) ;

    fParticleChange.ProposePolarization(State(fTransportEndSpin));

    G4double deltaTime = 0.0 ;

    G4double startTime = track.GetGlobalTime() ;

    if (State(fEndGlobalTimeComputed) == false) // means the track is leading the step
    {
        if(GetIT(track)->GetTrackingInfo()->IsLeadingStep() == false )
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "The track " << track.GetTrackID()
                                 << " should be the step leader";
            G4Exception("G4DNABrownianTransportation::AlongStepDoIt","G4DNABrownianTransportation002",
                        FatalErrorInArgument,exceptionDescription);
        }

        deltaTime = step.GetPostStepPoint()->GetGlobalTime() - startTime;

        fParticleChange.ProposeGlobalTime(  State(fCandidateEndGlobalTime) );
        fParticleChange.ProposePosition(State(fTransportEndPosition)) ;
    }
    else // the track should not be leading the step
    {
        deltaTime = State(fCandidateEndGlobalTime) - startTime ;

        if(GetIT(track)->GetTrackingInfo()->IsLeadingStep())
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "The track " << track.GetTrackID()
                                 << " should not be the step leader";
            G4Exception("G4DNABrownianTransportation::AlongStepDoIt","G4DNABrownianTransportation002",
                        FatalErrorInArgument,exceptionDescription);
        }

        fParticleChange.ProposePosition(State(fTransportEndPosition)) ;

        ///////////////
        fParticleChange.ProposeGlobalTime( State(fCandidateEndGlobalTime) ) ;
        ///////////////
    }

#ifdef G4VERBOSE
    //    DEBUG
    if(fVerboseLevel > 1)
    {
        G4cout<< GREEN_ON_BLUE<<"G4ITBrownianTransportation::AlongStepDoIt() :"
              << " trackID : " << track.GetTrackID()
              <<" Molecule name: "<< GetMolecule(track)-> GetName()
              << G4endl
              << "Current diffusion length : "<< G4BestUnit(step.GetStepLength(),"Length") <<" within time step : "
              << G4BestUnit((State(fCandidateEndGlobalTime) - step.GetPreStepPoint() -> GetGlobalTime()), "Time")
              << G4endl
              << "Propose diffusion length :" << (State(fTransportEndPosition) - track.GetPosition()).mag()
              << RESET
              << G4endl << G4endl;
    }
#endif

    G4double  restMass       = track.GetDynamicParticle()->GetMass() ;
    G4double deltaProperTime = deltaTime*( restMass/track.GetTotalEnergy() ) ;

    fParticleChange.ProposeProperTime(track.GetProperTime() + deltaProperTime) ;
    //fParticleChange. ProposeTrueStepLength( track.GetStepLength() ) ;

    // If the particle is caught looping or is stuck (in very difficult
    // boundaries) in a magnetic field (doing many steps)
    //   THEN this kills it ...
    //
    if ( State(fParticleIsLooping) )
    {
        G4double endEnergy= State(fTransportEndKineticEnergy);

        if( (endEnergy < fThreshold_Important_Energy)
                || (State(fNoLooperTrials) >= fThresholdTrials ) )
        {
            // Kill the looping particle
            //
            //            G4cout << "G4DNABrownianTransportation will killed the molecule"<< G4endl;
            fParticleChange.ProposeTrackStatus( fStopAndKill )  ;

            // 'Bare' statistics
            fSumEnergyKilled += endEnergy;
            if( endEnergy > fMaxEnergyKilled)
            {
                fMaxEnergyKilled= endEnergy;
            }

#ifdef G4VERBOSE
            if( (fVerboseLevel > 1) ||
                    ( endEnergy > fThreshold_Warning_Energy )  )
            {
                G4cout << " G4DNABrownianTransportation is killing track that is looping or stuck "
                       << G4endl
                       << "   This track has " << track.GetKineticEnergy() / MeV
                       << " MeV energy." << G4endl;
                G4cout << "   Number of trials = " << State(fNoLooperTrials)
                       << "   No of calls to AlongStepDoIt = " << noCalls
                       << G4endl;
            }
#endif
            State(fNoLooperTrials)=0;
        }
        else
        {
            State(fNoLooperTrials) ++;
#ifdef G4VERBOSE
            if( (fVerboseLevel > 2) )
            {
                G4cout << "   G4DNABrownianTransportation::AlongStepDoIt(): Particle looping -  "
                       << "   Number of trials = " << State(fNoLooperTrials)
                       << "   No of calls to  = " << noCalls
                       << G4endl;
            }
#endif
        }
    }
    else
    {
        State(fNoLooperTrials)=0;
    }

    // Another (sometimes better way) is to use a user-limit maximum Step size
    // to alleviate this problem ..

    return &fParticleChange ;
}
