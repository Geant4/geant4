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
// $Id: G4DNABrownianTransportation.cc 64374 2012-10-31 16:37:23Z gcosmo $
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

/// \brief { The transportation method implemented is the one from
///         Ermak-McCammon : J. Chem. Phys. 69, 1352 (1978)}

#include <CLHEP/Random/Stat.h>

#include "G4DNABrownianTransportation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Molecule.hh"
#include "G4RandomDirection.hh"
#include "G4ParticleTable.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4DNAMolecularMaterial.hh"

using namespace std;

#ifndef State
#define State(theXInfo) (fpBrownianState->theXInfo)
#endif

//#ifndef State
//#define State(theXInfo) (fTransportationState->theXInfo)
//#endif

//COLOR FOR DEBUGING
//#define RED_ON_WHITE  "\033[0;31m"
//#define GREEN "\033[32;40m"
#define GREEN_ON_BLUE "\033[1;32;44m"
#define RESET "\033[0m"

static double InvErf(double x)
{
    return CLHEP::HepStat::inverseErf(x);
}

static double InvErfc(double x)
{
    return CLHEP::HepStat::inverseErf(1.-x);
}

G4DNABrownianTransportation::G4DNABrownianTransportation(const G4String& aName, G4int verbosity) :
    G4ITTransportation(aName, verbosity),
    InitProcessState(fpBrownianState, fTransportationState)
{
    //ctor
    SetProcessSubType(61);
    verboseLevel = 1;
    fUseMaximumTimeBeforeReachingBoundary = true;
    fNistWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    fpWaterDensity  = 0;
}

G4DNABrownianTransportation::~G4DNABrownianTransportation()
{;}

G4DNABrownianTransportation::G4DNABrownianTransportation(const G4DNABrownianTransportation& right) :
    G4ITTransportation(right),
    InitProcessState(fpBrownianState, fTransportationState)
{
    //copy ctor
    SetProcessSubType(61);
    verboseLevel = right.verboseLevel;
    fUseMaximumTimeBeforeReachingBoundary = right.fUseMaximumTimeBeforeReachingBoundary;
    fNistWater = right.fNistWater;
    fpWaterDensity  = right.fpWaterDensity ;
}

G4DNABrownianTransportation& G4DNABrownianTransportation::operator=(const G4DNABrownianTransportation& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

G4DNABrownianTransportation::G4ITBrownianState::G4ITBrownianState() : G4ITTransportationState()
{
    fPathLengthWasCorrected = false;
}

void G4DNABrownianTransportation::StartTracking(G4Track* track)
{
    fpState = new G4ITBrownianState();
    SetInstantiateProcessState(false);
    G4ITTransportation::StartTracking(track);
}

void G4DNABrownianTransportation::BuildPhysicsTable(const G4ParticleDefinition& particle)
{
    if(verboseLevel > 0)
    {
        G4cout << G4endl << GetProcessName() << ":   for  "
               << setw(24) << particle.GetParticleName()
               << "\tSubType= " << GetProcessSubType()   << G4endl;
    }

    // Initialize water density pointer
    fpWaterDensity = G4DNAMolecularMaterial::Instance()->GetDensityTableFor(G4Material::GetMaterial("G4_WATER"));
}

void G4DNABrownianTransportation::ComputeStep(const G4Track& track,
                                              const G4Step& step,
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
        const G4VITProcess* ITProc = ((const G4VITProcess*) step.GetPostStepPoint()->GetProcessDefinedStep());
        bool makeException = true;

        if(ITProc && ITProc->ProposesTimeStep()) makeException = false;

        if(makeException)
        {

        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "ComputeStep is called while the track has the minimum interaction time";
        exceptionDescription << " so it should not recompute a timeStep ";
        G4Exception("G4DNABrownianTransportation::ComputeStep","G4DNABrownianTransportation001",
                    FatalErrorInArgument,exceptionDescription);
        }
    }

    State(fGeometryLimitedStep) = false;
    // TODO : generalize this process to all kind of brownian objects
    //    G4ITBrownianObject* ITBrown = GetITBrownianObject(track) ;
    //    G4double diffCoeff = ITBrown->GetDiffusionCoefficient(track.GetMaterial());
    G4Molecule* molecule = GetMolecule(track) ;
    G4double diffCoeff = molecule->GetDiffusionCoefficient();

    if(timeStep > 0)
    {
        spaceStep = DBL_MAX;

        while(spaceStep > State(endpointDistance))
            // Probably inefficient when the track is close close to boundaries
            // it goes with fUserMaximumTimeBeforeReachingBoundary == false
            // fUserMaximumTimeBeforeReachingBoundary == true, it should never loop
        {
            G4double x = G4RandGauss::shoot(0,sqrt(2*diffCoeff*timeStep));
            G4double y = G4RandGauss::shoot(0,sqrt(2*diffCoeff*timeStep));
            G4double z = G4RandGauss::shoot(0,sqrt(2*diffCoeff*timeStep));

            spaceStep = sqrt(x*x + y*y + z*z);
        }
        //        State(fTransportEndPosition).set(x + track.GetPosition().x(),
        //                                         y + track.GetPosition().y(),
        //                                         z + track.GetPosition().z());

        State(fTransportEndPosition)= spaceStep*step.GetPostStepPoint()->GetMomentumDirection() + track.GetPosition();
    }
    else
    {
        spaceStep = 0. ;
        State(fTransportEndPosition) =  track.GetPosition() ;
    }

    State(fCandidateEndGlobalTime) = step.GetPreStepPoint()->GetGlobalTime() + timeStep ;
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
//    if (material != fNistWater && material->GetBaseMaterial() != fNistWater)

   G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

   if(waterDensity == 0.0)
//  if (material == nistwater || material->GetBaseMaterial() == nistwater)
    {
        G4cout << "A track is outside water material : trackID"<< track.GetTrackID() << " (" << GetMolecule(track)->GetName()  <<")" << G4endl;
        G4cout << "Local Time : "<<  (track.GetLocalTime()) /s<<G4endl;
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
    State(fMomentumChanged) = true;
    fParticleChange.SetMomentumChanged(true) ;

    //    G4cout << "BROWN : Propose new direction :" << G4ThreeVector(xMomentum, yMomentum, zMomentum) << G4endl;

    // Alternative
    //fParticleChange.ProposeMomentumDirection(G4RandomDirection());

    return; // &fParticleChange is the final returned object
}


G4double G4DNABrownianTransportation::AlongStepGetPhysicalInteractionLength(
        const G4Track& track,
        G4double previousStepSize,
        G4double currentMinimumStep,
        G4double& currentSafety,
        G4GPILSelection* selection)
{
    G4double geometryStepLength = G4ITTransportation::AlongStepGetPhysicalInteractionLength(track,
                                                                                            previousStepSize,
                                                                                            currentMinimumStep,
                                                                                            currentSafety,
                                                                                            selection);

    G4double diffusionCoefficient = GetMolecule(track)->GetDiffusionCoefficient();
    //    G4double diffusionCoefficient = GetITBrownianObject(track)->GetDiffusionCoefficient(track.GetMaterial());

    if(State(fGeometryLimitedStep))
    {
        // 99 % of the space step distribution is lower than
        // d_99 = 8 * sqrt(D*t)
        // where t is the corresponding time step
        // so by inversion :
        if(fUseMaximumTimeBeforeReachingBoundary)
        {
            State(theInteractionTimeLeft)     = (geometryStepLength*geometryStepLength)/(64 * diffusionCoefficient);
        }
        else // Will use a random time
        {
            State(theInteractionTimeLeft) = 1/(4*diffusionCoefficient) * pow(geometryStepLength/InvErfc(G4UniformRand()),2);
        }

        State(fCandidateEndGlobalTime) = track.GetGlobalTime() + State(theInteractionTimeLeft);
        State(fPathLengthWasCorrected) = false;
    }
    else
    {
        geometryStepLength  = 2*sqrt(diffusionCoefficient*State(theInteractionTimeLeft) ) *InvErf(G4UniformRand());
        State(fPathLengthWasCorrected) = true;
        State(endpointDistance) = geometryStepLength;
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
    G4ITTransportation::AlongStepDoIt(track,step);
    Diffusion(track);
    return &fParticleChange;
}
