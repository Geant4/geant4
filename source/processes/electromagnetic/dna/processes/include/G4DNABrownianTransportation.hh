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
// $Id: G4DNABrownianTransportation.hh 90232 2015-05-21 08:54:54Z gcosmo $
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4ITBROWNIANTRANSPORTATION_H
#define G4ITBROWNIANTRANSPORTATION_H

#include "G4ITTransportation.hh"

class G4SafetyHelper;
class G4Molecule;

// experimental
class G4BrownianAction
{
public:
  G4BrownianAction(){;}
  virtual ~G4BrownianAction(){;}

//  virtual G4double GetDiffusionCoefficient(G4Material*,
//                                           G4Molecule*) { return 0;}

  // If returns true: track is killed
  virtual void Transport(const G4Track&,
                         G4ParticleChangeForTransport&) = 0;
};


/* \brief {The transportation method implemented is the one from
 *         Ermak-McCammon : J. Chem. Phys. 69, 1352 (1978).
 *         To compute time and space intervals to reach a volume boundary,
 *         there are two alternative methods proposed by this process.
 *
 *         ** Method 1 selects a minimum distance to the next
 *         boundary using to the following formula:
 *
 *         --> t_min = (safety* safety) / (8 * diffusionCoefficient);
 *         this corresponds to 5% probability of the Brownian particle to cross
 *         the boundary - isotropic distance to nearest boundary (safety) is used
 *
 *         OR if the flag "speed me up" is on:
 *
 *         --> t_min = (geometryStepLength * geometryStepLength) / diffusionCoefficient;
 *         this corresponds to 50% probability of the Brownian particle to cross
 *         the boundary - distance along current direction to nearest boundary is used
 *
 *         NB: By default, method 1 with the flag "speed me up is used".
 *         In addition, one may want to used the minimum time step limit defined
 *         in G4Scheduler through the G4UserTimeStepAction. If so, speed level might
 *         be set to 2. But minimum time steps have to be set in the user class.
 *
 *         ** Method 2 can randomly compute the time to the next boundary using the
 *         following formula:
 *
 *         t_random = 1 / (4 * diffusionCoefficient)* pow(geometryStepLength /
 *                         InvErfc(G4UniformRand()),2);
 *         For release 10.1, this is using the 1D cumulative density function.
 *
 *         At each diffusion step, the direction of the particle is randomly selected.
 *         For now, the geometryStepLength corresponds to the distance to the
 *         nearest boundary along the direction of diffusion which selected randomly.
 *
 *         Method 2 is currently deactivated by default.
 *         }
 */

class G4DNABrownianTransportation : public G4ITTransportation
{
public:
  G4DNABrownianTransportation(const G4String& aName =
                              "DNABrownianTransportation",
                              G4int verbosityLevel = 0);
  G4IT_ADD_CLONE(G4VITProcess,G4DNABrownianTransportation)
  virtual ~G4DNABrownianTransportation();
  G4DNABrownianTransportation(const G4DNABrownianTransportation& other);
  G4DNABrownianTransportation& operator=(const G4DNABrownianTransportation& other);

  inline void SetBrownianAction(G4BrownianAction*);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual void StartTracking(G4Track* aTrack);

  virtual void ComputeStep(const G4Track&,
                           const G4Step&,
                           const double,
                           double&);

  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track& /*track*/,
                                                         G4double /*previousStepSize*/,
                                                         G4double /*currentMinimumStep*/,
                                                         G4double& /*currentSafety*/,
                                                         G4GPILSelection* /*selection*/);
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step&);

  virtual G4VParticleChange* AlongStepDoIt(const G4Track& track, const G4Step&);

  // Boundary is crossed at time at which:
  // * either 5% of the distribution might be over boundary - the last position
  //   is adjusted on boundary
  // * or if speedUp (from level 1) is activated - 50% of the distribution might
  //   be over boundary, the particles are also allowed to jump over boundary
  inline void UseMaximumTimeBeforeReachingBoundary(bool flag = true)
  {
    fUseMaximumTimeBeforeReachingBoundary = flag;
  }

  // Random sampling time at which boundary is crossed
  // WARNING: For release 10.1, this is a 1D approximation for sampling time
  // but 3D for diffusion
  // If speed up IS activated, particles are allowed jump over barrier
  inline void UseCumulativeDensitFunction(bool flag = true)
  {
    if(flag == true)
    {
      fUseMaximumTimeBeforeReachingBoundary = false;
      return;
    }
    fUseMaximumTimeBeforeReachingBoundary = true;
  }

  // Use limiting time steps defined in the scheduler
  inline void UseLimitingTimeSteps(bool flag = true)
  {
    fUseSchedulerMinTimeSteps = flag;
  }

  inline void SpeedLevel(int level)
  {
    if(level < 0) level =0;
    else if(level > 2) level = 2;

    switch(level)
    {
      case 0:
        fSpeedMeUp = false;
        fUseSchedulerMinTimeSteps = false;
        return;

      case 1:
        fSpeedMeUp = true;
        fUseSchedulerMinTimeSteps = false;
        return;

      case 2:
        //======================================================================
        // NB: BE AWARE THAT IF NO MIN TIME STEPS NO TIME STEPS HAVE BEEN
        // PROVIDED TO G4Scheduler THIS LEVEL MIGHT BE SLOWER THAN LEVEL 1
        //======================================================================
        fSpeedMeUp = true;
        fUseSchedulerMinTimeSteps = true;
        return;
    }
  }

protected:

  G4double ComputeGeomLimit(const G4Track& track,
                            G4double& presafety,
                            G4double limit);

  void Diffusion(const G4Track& track);

  //________________________________________________________________
  // Process information
  struct G4ITBrownianState : public G4ITTransportationState
  {
  public:
    G4ITBrownianState();
    virtual ~G4ITBrownianState()
    {
      ;
    }
    virtual G4String GetType()
    {
      return "G4ITBrownianState";
    }

    G4bool fPathLengthWasCorrected;
    G4bool fTimeStepReachedLimit;
    G4bool fComputeLastPosition;
    G4double  fRandomNumber;
  };

  G4bool fUseMaximumTimeBeforeReachingBoundary;
  G4Material* fNistWater;

  G4bool fUseSchedulerMinTimeSteps;
  G4double  fInternalMinTimeStep;
  G4bool fSpeedMeUp;

  // Water density table
  const std::vector<G4double>* fpWaterDensity;

  G4BrownianAction* fpBrownianAction;
};


inline void G4DNABrownianTransportation::SetBrownianAction(G4BrownianAction* brownianAction)
{
  fpBrownianAction = brownianAction;
}

#endif // G4ITBROWNIANTRANSPORTATION_H
