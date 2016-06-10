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
// $Id: G4DNABrownianTransportation.cc 94218 2015-11-09 08:24:48Z gcosmo $
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

#include <G4Scheduler.hh>
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
#include "G4ITNavigator.hh"
#include "G4ITSafetyHelper.hh" // Not used yet
#include "G4TrackingInformation.hh"

using namespace std;

#ifndef State
#define State(theXInfo) (GetState<G4ITBrownianState>()->theXInfo)
#endif

//#ifndef State
//#define State(theXInfo) (fTransportationState->theXInfo)
//#endif

//#define USE_COLOR 1

#ifdef USE_COLOR
#define RED  "\033[0;31m"
#define LIGHT_RED  "\33[1;31m"
#define GREEN "\033[32;40m"
#define GREEN_ON_BLUE "\033[1;32;44m"
#define RESET_COLOR "\033[0m"
#else
#define RED  ""
#define LIGHT_RED  ""
#define GREEN ""
#define GREEN_ON_BLUE ""
#define RESET_COLOR ""
#endif

//#define DEBUG_MEM 1

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
using G4MemStat::MemStat;
#endif

static double InvErf(double x)
{
  return CLHEP::HepStat::inverseErf(x);
}

static double InvErfc(double x)
{
  return CLHEP::HepStat::inverseErf(1. - x);
}

//static double Erf(double x)
//{
//  return CLHEP::HepStat::erf(x);
//}

static double Erfc(double x)
{
  return 1 - CLHEP::HepStat::erf(1. - x);
}

#ifndef State
#define State(theXInfo) (GetState<G4ITTransportationState>()->theXInfo)
#endif

G4DNABrownianTransportation::G4DNABrownianTransportation(const G4String& aName,
                                                         G4int verbosity) :
    G4ITTransportation(aName, verbosity)
{

  fVerboseLevel = 0;

  fpState.reset(new G4ITBrownianState());

  //ctor
  SetProcessSubType(61);

  fNistWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  fpWaterDensity = 0;

  fUseMaximumTimeBeforeReachingBoundary = true;
  fUseSchedulerMinTimeSteps = false;
  fSpeedMeUp = true;

  fInternalMinTimeStep = 1*picosecond;
  fpBrownianAction = 0;
}

G4DNABrownianTransportation::~G4DNABrownianTransportation()
{
  if(fpBrownianAction) delete fpBrownianAction;
}

G4DNABrownianTransportation::G4DNABrownianTransportation(const G4DNABrownianTransportation& right) :
    G4ITTransportation(right)
{
  //copy ctor
  SetProcessSubType(61);
  fUseMaximumTimeBeforeReachingBoundary = right
      .fUseMaximumTimeBeforeReachingBoundary;
  fUseSchedulerMinTimeSteps = right.fUseSchedulerMinTimeSteps;
  fNistWater = right.fNistWater;
  fpWaterDensity = right.fpWaterDensity;
  fInternalMinTimeStep = right.fInternalMinTimeStep;
  fSpeedMeUp = right.fSpeedMeUp;
  fpBrownianAction = right.fpBrownianAction;
}

G4DNABrownianTransportation& G4DNABrownianTransportation::operator=(const G4DNABrownianTransportation& rhs)
{
  if(this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

G4DNABrownianTransportation::G4ITBrownianState::G4ITBrownianState() :
    G4ITTransportationState()
{
  fPathLengthWasCorrected = false;
  fTimeStepReachedLimit = false;
  fComputeLastPosition = false;
  fRandomNumber = -1;
}

void G4DNABrownianTransportation::StartTracking(G4Track* track)
{
  fpState.reset(new G4ITBrownianState());
//	G4cout << "G4DNABrownianTransportation::StartTracking : "
//  "Initialised track State" << G4endl;
  SetInstantiateProcessState(false);
  G4ITTransportation::StartTracking(track);
}

void G4DNABrownianTransportation::BuildPhysicsTable(const G4ParticleDefinition& particle)
{
  if(verboseLevel > 0)
  {
    G4cout << G4endl<< GetProcessName() << ":   for  "
    << setw(24) << particle.GetParticleName()
    << "\tSubType= " << GetProcessSubType() << G4endl;
  }
  // Initialize water density pointer
  fpWaterDensity = G4DNAMolecularMaterial::Instance()->
  GetDensityTableFor(G4Material::GetMaterial("G4_WATER"));

  fpSafetyHelper->InitialiseHelper();
  G4ITTransportation::BuildPhysicsTable(particle);
}

void G4DNABrownianTransportation::ComputeStep(const G4Track& track,
                                              const G4Step& step,
                                              const double timeStep,
                                              double& spaceStep)
{
  // G4cout << "G4ITBrownianTransportation::ComputeStep" << G4endl;

  /* If this method is called, this step
   * cannot be geometry limited.
   * In case the step is limited by the geometry,
   * this method should not be called.
   * The fTransportEndPosition calculated in
   * the method AlongStepIL should be taken
   * into account.
   * In order to do so, the flag IsLeadingStep
   * is on. Meaning : this track has the minimum
   * interaction length over all others.
   */
  if(GetIT(track)->GetTrackingInfo()->IsLeadingStep())
  {
    const G4VITProcess* ITProc = ((const G4VITProcess*) step.GetPostStepPoint()
        ->GetProcessDefinedStep());
    bool makeException = true;

    if(ITProc && ITProc->ProposesTimeStep()) makeException = false;

    if(makeException)
    {

      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "ComputeStep is called while the track has"
                              "the minimum interaction time";
      exceptionDescription << " so it should not recompute a timeStep ";
      G4Exception("G4DNABrownianTransportation::ComputeStep",
                  "G4DNABrownianTransportation001", FatalErrorInArgument,
                  exceptionDescription);
    }
  }

  State(fGeometryLimitedStep) = false;

  G4Molecule* molecule = GetMolecule(track);

  if(timeStep > 0)
  {
    spaceStep = DBL_MAX;

    // TODO : generalize this process to all kind of Brownian objects
    G4double diffCoeff = molecule->GetDiffusionCoefficient(track.GetMaterial(),
                                                           track.GetMaterial()->GetTemperature());

    static double sqrt_2 = sqrt(2.);
    double sqrt_Dt = sqrt(diffCoeff*timeStep);
    double sqrt_2Dt = sqrt_2*sqrt_Dt;

    if(State(fTimeStepReachedLimit)== true)
    {
      //========================================================================
      State(fGeometryLimitedStep) = true;// important
      //========================================================================
      spaceStep = State(fEndPointDistance);
   //   G4cout << "State(fTimeStepReachedLimit)== true" << G4endl;
    }
    else
    {
      G4double x = G4RandGauss::shoot(0,sqrt_2Dt);
      G4double y = G4RandGauss::shoot(0,sqrt_2Dt);
      G4double z = G4RandGauss::shoot(0,sqrt_2Dt);

      spaceStep = sqrt(x*x + y*y + z*z);

      if(spaceStep >= State(fEndPointDistance))
      {
        //G4cout << "spaceStep >= State(fEndPointDistance)" << G4endl;
        //======================================================================
        State(fGeometryLimitedStep) = true;// important
        //======================================================================
/*
        if(fSpeedMeUp)
        {
          G4cout << "SpeedMeUp" << G4endl;
        }
        else
*/
        if(fUseSchedulerMinTimeSteps == false)// jump over barrier NOT used
        {
#ifdef G4VERBOSE
          if (fVerboseLevel > 1)
          {
            G4cout << GREEN_ON_BLUE
            << "G4ITBrownianTransportation::ComputeStep() : "
            << "Step was limited to boundary"
            << RESET_COLOR
            << G4endl;
          }
#endif
          //TODO
          if(State(fRandomNumber)>=0) // CDF is used
          {
            /*
             //=================================================================
             State(fGeometryLimitedStep) = true;// important
             //=================================================================
             spaceStep = State(fEndPointDistance);
             */

            //==================================================================
            // BE AWARE THAT THE TECHNIQUE USED BELOW IS A 1D APPROXIMATION
            // Cumulative density function for the 3D case is not yet
            // implemented
            //==================================================================
//            G4cout << GREEN_ON_BLUE
//                   << "G4ITBrownianTransportation::ComputeStep() : "
//                   << "A random number was selected"
//                   << RESET_COLOR
//                   << G4endl;
            double value = State(fRandomNumber)+(1-State(fRandomNumber))*G4UniformRand();
            double invErfc = InvErfc(value);
            spaceStep = invErfc*2*sqrt_Dt;

            if(State(fTimeStepReachedLimit)== false)
            {
              //================================================================
              State(fGeometryLimitedStep) = false;// important
              //================================================================
            }
            //==================================================================
            // DEBUG
//            if(spaceStep > State(fEndPointDistance))
//            {
//              G4cout << "value = " << value << G4endl;
//              G4cout << "invErfc = " << invErfc << G4endl;
//              G4cout << "spaceStep = " << G4BestUnit(spaceStep, "Length")
//              << G4endl;
//              G4cout << "end point distance= " << G4BestUnit(State(fEndPointDistance), "Length")
//              << G4endl;
//            }
//
//            assert(spaceStep <= State(fEndPointDistance));
            //==================================================================

          }
          else if(fUseMaximumTimeBeforeReachingBoundary == false) // CDF is used
          {
            double min_randomNumber = Erfc(State(fEndPointDistance)/2*sqrt_Dt);
            double value = min_randomNumber+(1-min_randomNumber)*G4UniformRand();
            double invErfc = InvErfc(value);
            spaceStep = invErfc*2*sqrt_Dt;
            if(spaceStep >= State(fEndPointDistance))
            {
              //================================================================
              State(fGeometryLimitedStep) = true;// important
              //================================================================
            }
            else if(State(fTimeStepReachedLimit)== false)
            {
              //================================================================
              State(fGeometryLimitedStep) = false;// important
              //================================================================
            }
          }
          else // CDF is NOT used
          {
            //==================================================================
            State(fGeometryLimitedStep) = true;// important
            //==================================================================
            spaceStep = State(fEndPointDistance);
            //TODO

            /*
            //==================================================================
            // 1D approximation to place the brownian between its starting point
            // and the geometry boundary
            //==================================================================
            double min_randomNumber = Erfc(State(fEndPointDistance)/2*sqrt_Dt);
            double value = State(fRandomNumber)+(1-State(fRandomNumber))*G4UniformRand();
            double invErfc = InvErfc(value*G4UniformRand());
            spaceStep = invErfc*2*sqrt_Dt;
            State(fGeometryLimitedStep) = false;
            */
          }
        }

        State(fTransportEndPosition)= spaceStep*
//             step.GetPostStepPoint()->GetMomentumDirection() 
             track.GetMomentumDirection()
             + track.GetPosition();
      }
      else
      {
        //======================================================================
        State(fGeometryLimitedStep) = false;// important
        //======================================================================
        State(fTransportEndPosition)= spaceStep*step.GetPostStepPoint()->
        GetMomentumDirection() + track.GetPosition();
      }
    }
  }
  else
  {
    spaceStep = 0.;
    State(fTransportEndPosition) = track.GetPosition();
    State(fGeometryLimitedStep) = false;
  }

  State(fCandidateEndGlobalTime) = step.GetPreStepPoint()->GetGlobalTime()
      + timeStep;
  State(fEndGlobalTimeComputed) = true;

#ifdef G4VERBOSE
  //    DEBUG
  if (fVerboseLevel > 1)
  {
    G4cout << GREEN_ON_BLUE
           << "G4ITBrownianTransportation::ComputeStep() : "
           << " trackID : " << track.GetTrackID() << " : Molecule name: "
           << molecule->GetName() << G4endl
           << "Initial position:" << G4BestUnit(track.GetPosition(), "Length")
           << G4endl
           << "Initial direction:" << track.GetMomentumDirection() << G4endl
           << "Final position:" << G4BestUnit(State(fTransportEndPosition), "Length")
           << G4endl
           << "Initial magnitude:" << G4BestUnit(track.GetPosition().mag(), "Length")
           << G4endl
           << "Final magnitude:" << G4BestUnit(State(fTransportEndPosition).mag(), "Length")
           << G4endl
           << "Diffusion length : "
           << G4BestUnit(spaceStep, "Length")
           << " within time step : " << G4BestUnit(timeStep,"Time")
           << G4endl
           << "State(fTimeStepReachedLimit)= " << State(fTimeStepReachedLimit) << G4endl
           << "State(fGeometryLimitedStep)=" << State(fGeometryLimitedStep) << G4endl
           << "End point distance was: " << G4BestUnit(State(fEndPointDistance), "Length")
           << G4endl
           << RESET_COLOR
           << G4endl<< G4endl;
  }
#endif

//==============================================================================
// DEBUG
//assert(spaceStep <  State(fEndPointDistance)
//  || (spaceStep >= State(fEndPointDistance) && State(fGeometryLimitedStep)));
//assert(track.GetMomentumDirection() == State(fTransportEndMomentumDir));
//==============================================================================
}

G4VParticleChange* G4DNABrownianTransportation::PostStepDoIt(const G4Track& track,
                                                             const G4Step& step)
{
  G4ITTransportation::PostStepDoIt(track, step);

#ifdef G4VERBOSE
  //    DEBUG
  if (fVerboseLevel > 1)
  {
    G4cout << GREEN_ON_BLUE << "G4ITBrownianTransportation::PostStepDoIt() :"
           << " trackID : " << track.GetTrackID() << " Molecule name: "
           << GetMolecule(track)->GetName() << G4endl;
    G4cout << "Diffusion length : "
           << G4BestUnit(step.GetStepLength(), "Length")
           <<" within time step : " << G4BestUnit(step.GetDeltaTime(),"Time")
           << "\t Current global time : "
           << G4BestUnit(track.GetGlobalTime(),"Time")
           << RESET_COLOR
           << G4endl<< G4endl;
  }
#endif
  return &fParticleChange;
}

void G4DNABrownianTransportation::Diffusion(const G4Track& track)
{

#ifdef DEBUG_MEM
  MemStat mem_first = MemoryUsage();
#endif

#ifdef G4VERBOSE
  // DEBUG
  if (fVerboseLevel > 1)
  {
    G4cout << GREEN_ON_BLUE << setw(18)
           << "G4DNABrownianTransportation::Diffusion :" << setw(8)
           << GetIT(track)->GetName() << "\t trackID:" << track.GetTrackID()
           << "\t" << " Global Time = "
           << G4BestUnit(track.GetGlobalTime(), "Time")
           << RESET_COLOR
           << G4endl
           << G4endl;
  }
#endif

/*
  fParticleChange.ProposePosition(State(fTransportEndPosition));
  //fParticleChange.ProposeEnergy(State(fTransportEndKineticEnergy));
  fParticleChange.SetMomentumChanged(State(fMomentumChanged));

  fParticleChange.ProposeGlobalTime(State(fCandidateEndGlobalTime));
  fParticleChange.ProposeLocalTime(State(fCandidateEndGlobalTime));
  fParticleChange.ProposeTrueStepLength(track.GetStepLength());
*/
  G4Material* material = track.GetMaterial();

  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

  if (waterDensity == 0.0)
  {
    if(fpBrownianAction)
    {
      // Let the user Brownian action class decide what to do
      fpBrownianAction->Transport(track,
                                  fParticleChange);
      return;
    }
    else
    {
#ifdef G4VERBOSE
      if(fVerboseLevel)
      {
        G4cout << "A track is outside water material : trackID = "
        << track.GetTrackID() << " (" << GetMolecule(track)->GetName() <<")"
        << G4endl;
        G4cout << "Local Time : " << G4BestUnit(track.GetGlobalTime(), "Time")
        << G4endl;
        G4cout << "Step Number :" << track.GetCurrentStepNumber() << G4endl;
      }
#endif
      fParticleChange.ProposeEnergy(0.);
      fParticleChange.ProposeTrackStatus(fStopAndKill);
      return;// &fParticleChange is the final returned object
    }
  }


   #ifdef DEBUG_MEM
   MemStat mem_intermediaire = MemoryUsage();
   mem_diff = mem_intermediaire-mem_first;
   G4cout << "\t\t\t >> || MEM || In G4DNABrownianTransportation::Diffusion "
   "after dealing with waterDensity for "<< track.GetTrackID()
   << ", diff is : " << mem_diff << G4endl;
   #endif

  fParticleChange.ProposeMomentumDirection(G4RandomDirection());
  State(fMomentumChanged) = true;
  fParticleChange.SetMomentumChanged(true);
   //
   #ifdef DEBUG_MEM
   mem_intermediaire = MemoryUsage();
   mem_diff = mem_intermediaire-mem_first;
   G4cout << "\t\t\t >> || MEM || In G4DNABrownianTransportation::"
          "After proposing new direction to fParticleChange for "
          << track.GetTrackID() << ", diff is : " << mem_diff << G4endl;
   #endif

  return;// &fParticleChange is the final returned object
}

// NOT USED
G4double G4DNABrownianTransportation::ComputeGeomLimit(const G4Track& track,
                                                       G4double& presafety,
                                                       G4double limit)
{
  G4double res = DBL_MAX;
  if(track.GetVolume() != fpSafetyHelper->GetWorldVolume())
  {
    G4TrackStateManager& trackStateMan = GetIT(track)->GetTrackingInfo()
    ->GetTrackStateManager();
    fpSafetyHelper->LoadTrackState(trackStateMan);
    res = fpSafetyHelper->CheckNextStep(
        track.GetStep()->GetPreStepPoint()->GetPosition(),
        track.GetMomentumDirection(),
        limit, presafety);
    fpSafetyHelper->ResetTrackState();
  }
  return res;
}

G4double G4DNABrownianTransportation::AlongStepGetPhysicalInteractionLength(const G4Track& track,
                                                                            G4double previousStepSize,
                                                                            G4double currentMinimumStep,
                                                                            G4double& currentSafety,
                                                                            G4GPILSelection* selection)
{
#ifdef G4VERBOSE
  if(fVerboseLevel)
  {
    G4cout << " G4DNABrownianTransportation::AlongStepGetPhysicalInteractionLength - track ID: "
           << track.GetTrackID() << G4endl;
    G4cout << "In volume : " << track.GetVolume()->GetName()
           << " position : " << G4BestUnit(track.GetPosition() , "Length") << G4endl;
  }
#endif

  G4double geometryStepLength =
      G4ITTransportation::AlongStepGetPhysicalInteractionLength(
          track, previousStepSize, currentMinimumStep, currentSafety,
          selection);

  if(geometryStepLength==0)
  { 
//    G4cout << "geometryStepLength==0" << G4endl;
    if(State(fGeometryLimitedStep))
    {
//      G4cout << "if(State(fGeometryLimitedStep))" << G4endl;
      G4TouchableHandle newTouchable = new G4TouchableHistory;

      newTouchable->UpdateYourself(State(fCurrentTouchableHandle)->GetVolume(),
                                   State(fCurrentTouchableHandle)->GetHistory());

       fLinearNavigator->SetGeometricallyLimitedStep();
       fLinearNavigator->LocateGlobalPointAndUpdateTouchableHandle(
           track.GetPosition(), track.GetMomentumDirection(),
           newTouchable, true);

       if(newTouchable->GetVolume() == 0)
       {
          return 0;
       }

       State(fCurrentTouchableHandle) = newTouchable;

       //=======================================================================
       // TODO: speed up navigation update
//       geometryStepLength = fLinearNavigator->ComputeStep(track.GetPosition(),
//                                     track.GetMomentumDirection(),
//                                     currentMinimumStep,
//                                     currentSafety);
       //=======================================================================


       //=======================================
       // Longer but safer ...
       geometryStepLength =
             G4ITTransportation::AlongStepGetPhysicalInteractionLength(
                 track, previousStepSize, currentMinimumStep, currentSafety,
                 selection);

      }
  }

  //============================================================================
  // DEBUG
 // G4cout << "geometryStepLength: " << G4BestUnit(geometryStepLength, "Length")
 //        << " | trackID: " << track.GetTrackID()
 //        << G4endl;
  //============================================================================

  G4double diffusionCoefficient = 0;

  /*
  G4Material* material = track.GetMaterial();
  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

  if (waterDensity == 0.0)
  {
    if(fpBrownianAction)
    {
      diffusionCoefficient = fpBrownianAction->GetDiffusionCoefficient(material,
                                                                       GetMolecule(track));
    }

    if(diffusionCoefficient <= 0)
    {
      State(fGeometryLimitedStep) = false;
      State(theInteractionTimeLeft) = 0;
      State(fTransportEndPosition) = track.GetPosition();
      return 0;
    }

  }
  else
    */
  //{
    diffusionCoefficient = GetMolecule(track)->GetDiffusionCoefficient();
  //}

  State(fComputeLastPosition) = false;
  State(fTimeStepReachedLimit) = false;

  if (State(fGeometryLimitedStep))
  {
    // 95 % of the space step distribution is lower than
    // d_95 = 2 * sqrt(2*D*t)
    // where t is the corresponding time step
    // so by inversion :
    if (fUseMaximumTimeBeforeReachingBoundary)
    {
      if(fSpeedMeUp)
      {
      State(theInteractionTimeLeft) = (geometryStepLength * geometryStepLength)
          / (diffusionCoefficient); // d_50 - use straight line
      }
      else
      {
         State(theInteractionTimeLeft) = (currentSafety * currentSafety)
                  / (diffusionCoefficient); // d_50 - use safety

         //=====================================================================
         // State(theInteractionTimeLeft) = (currentSafety * currentSafety)
         //          / (8 * diffusionCoefficient); // d_95
         //=====================================================================
      }
      State(fComputeLastPosition) = true;
    }
    else
    // Will use a random time - this is precise but long to compute in certain
    // circumstances (many particles - small volumes)
    {
      State(fRandomNumber) = G4UniformRand();
      State(theInteractionTimeLeft) = 1 / (4 * diffusionCoefficient)
          * pow(geometryStepLength / InvErfc(State(fRandomNumber)),2);

      State(fTransportEndPosition) = geometryStepLength*
          track.GetMomentumDirection() + track.GetPosition();
    }

    if (fUseSchedulerMinTimeSteps)
    {
      double minTimeStepAllowed = G4VScheduler::Instance()->GetLimitingTimeStep();
      //========================================================================
      // TODO
//      double currentMinTimeStep = G4VScheduler::Instance()->GetTimeStep();
      //========================================================================

      if (State(theInteractionTimeLeft) < minTimeStepAllowed)
      {
        State(theInteractionTimeLeft) = minTimeStepAllowed;
        State(fTimeStepReachedLimit) = true;
        State(fComputeLastPosition) = true;
      }
    }
    else if(State(theInteractionTimeLeft) < fInternalMinTimeStep)
      // TODO: find a better way when fForceLimitOnMinTimeSteps is not used
    {
      State(fTimeStepReachedLimit) = true;
      State(theInteractionTimeLeft) = fInternalMinTimeStep;
      if (fUseMaximumTimeBeforeReachingBoundary)
      {
        State(fComputeLastPosition) = true;
      }
    }

    State(fCandidateEndGlobalTime) =
        track.GetGlobalTime() + State(theInteractionTimeLeft);

    State(fEndGlobalTimeComputed) = true; // MK: ADDED ON 20/11/2014

    State(fPathLengthWasCorrected) = false;
  }
  else
  {
    // Transform geometrical step
    geometryStepLength = 2
        * sqrt(diffusionCoefficient * State(theInteractionTimeLeft))
        * InvErf(G4UniformRand());
    State(fPathLengthWasCorrected) = true;
    //State(fEndPointDistance) = geometryStepLength;
    State(fTransportEndPosition) = geometryStepLength*
              track.GetMomentumDirection() + track.GetPosition();
  }

#ifdef G4VERBOSE
  //    DEBUG
  if (fVerboseLevel > 1)
  {
  G4cout << GREEN_ON_BLUE
         << "G4DNABrownianTransportation::AlongStepGetPhysicalInteractionLength = "
         << G4BestUnit(geometryStepLength, "Length")
         << " | trackID = "
         << track.GetTrackID()
         << RESET_COLOR
         << G4endl;
  }
#endif

// assert(geometryStepLength <  State(fEndPointDistance)
//  || (geometryStepLength >= State(fEndPointDistance) && State(fGeometryLimitedStep)));

  return geometryStepLength;
}

//////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)
G4VParticleChange*
G4DNABrownianTransportation::AlongStepDoIt(const G4Track& track,
                                           const G4Step& step)
{
#ifdef DEBUG_MEM
  MemStat mem_first, mem_second, mem_diff;
#endif

#ifdef DEBUG_MEM
  mem_first = MemoryUsage();
#endif

  if (GetIT(track)->GetTrackingInfo()->IsLeadingStep()
      && State(fComputeLastPosition))
  {
    //==========================================================================
    // DEBUG
    //
//     assert(fabs(State(theInteractionTimeLeft)-
//                G4VScheduler::Instance()->GetTimeStep()) < DBL_EPSILON);
    //==========================================================================

    double spaceStep = DBL_MAX;

    if(State(theInteractionTimeLeft) <= fInternalMinTimeStep)
    {
      spaceStep = State(fEndPointDistance);
      State(fGeometryLimitedStep) = true;
    }
    else
    {
      G4double diffusionCoefficient = GetMolecule(track)->
      GetDiffusionCoefficient();

      double sqrt_2Dt= sqrt(2 * diffusionCoefficient * State(theInteractionTimeLeft));
      G4double x = G4RandGauss::shoot(0, sqrt_2Dt);
      G4double y = G4RandGauss::shoot(0, sqrt_2Dt);
      G4double z = G4RandGauss::shoot(0, sqrt_2Dt);

      spaceStep = sqrt(x * x + y * y + z * z);

      if(spaceStep >= State(fEndPointDistance))
      {
        State(fGeometryLimitedStep) = true;
       if (
           //fSpeedMeUp == false&&
           fUseSchedulerMinTimeSteps == false
           && spaceStep >= State(fEndPointDistance))
       {
         spaceStep = State(fEndPointDistance);
       }
      }
      else
      {
        State(fGeometryLimitedStep) = false;
      }
    }

//    assert( (spaceStep <  State(fEndPointDistance) && State(fGeometryLimitedStep) == false)
//|| (spaceStep >= State(fEndPointDistance) && State(fGeometryLimitedStep)));

    // Calculate final position
    //
    State(fTransportEndPosition) = track.GetPosition()
               + spaceStep * track.GetMomentumDirection();
  }

  if(fVerboseLevel)
  {
    G4cout << GREEN_ON_BLUE
            << "G4DNABrownianTransportation::AlongStepDoIt: "
                "GeometryLimitedStep = "
            << State(fGeometryLimitedStep)
            << RESET_COLOR
            << G4endl;
  }

//  static G4ThreadLocal G4int noCalls = 0;
//  noCalls++;

//  fParticleChange.Initialize(track);

  G4ITTransportation::AlongStepDoIt(track, step);

#ifdef DEBUG_MEM
  MemStat mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire-mem_first;
  G4cout << "\t\t\t >> || MEM || After calling G4ITTransportation::"
  "AlongStepDoIt for "<< track.GetTrackID() << ", diff is : "
  << mem_diff << G4endl;
#endif

  if(track.GetStepLength() != 0 // && State(fGeometryLimitedStep)
      //========================================================================
      // TODO: avoid changing direction after too small time steps
//    && (G4VScheduler::Instance()->GetTimeStep() > fInternalMinTimeStep
//        || fSpeedMeUp == false)
      //========================================================================
        )
  {
    Diffusion(track);
  }
  //else
  //{
  //  fParticleChange.ProposeMomentumDirection(State(fTransportEndMomentumDir));
  //}
/*
  if (State(fParticleIsLooping))
  {
    if ((State(fNoLooperTrials)>= fThresholdTrials))
    {
      fParticleChange.ProposeTrackStatus(fStopAndKill);
      State(fNoLooperTrials) = 0;
#ifdef G4VERBOSE
      if ((fVerboseLevel > 1))
      {
        G4cout
            << " G4DNABrownianTransportation is killing track that is looping or stuck "
            << G4endl;
        G4cout << "   Number of trials = " << State(fNoLooperTrials)
        << "   No of calls to AlongStepDoIt = " << noCalls
        << G4endl;
      }
#endif
    }
    else
    {
      State(fNoLooperTrials)++;
    }
  }
  else
  {
    State(fNoLooperTrials)=0;
  }
*/
#ifdef DEBUG_MEM
  mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire-mem_first;
  G4cout << "\t\t\t >> || MEM || After calling G4DNABrownianTransportation::"
  "Diffusion for "<< track.GetTrackID() << ", diff is : "
  << mem_diff << G4endl;
#endif

  return &fParticleChange;
}

