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
// $Id: G4ITModelProcessor.hh 100802 2016-11-02 14:55:27Z gcosmo $
//
// Author: Mathieu Karamitros

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

#ifndef G4ITMODELPROCESSOR_H
#define G4ITMODELPROCESSOR_H

#include <vector>
#include "G4ITReactionChange.hh"
#include "G4ITType.hh"
#include "G4ITModelHandler.hh"
#include "G4ITStepStatus.hh"

class G4VITTimeStepComputer;
class G4VITReactionProcess;
class G4ITModelHandler;
class G4ITReactionSet;
class G4UserTimeStepAction;
class G4ITTrackingManager;
class G4ITTrackHolder;

//#ifndef compTrackPerID__
//#define compTrackPerID__
//  struct compTrackPerID
//  {
//    bool operator()(G4Track* rhs, G4Track* lhs) const
//    {
//      return rhs->GetTrackID() < lhs->GetTrackID();
//    }
//  };
//#endif

/**
 * The G4ITModelProcessor will call the two processes defined in G4VITModel.
 * This processes act at the beginning and end of each step.
 * The first one, the TimeStepper will calculate a time step to propagate all
 * the track and eventually it can return some tracks that can likely react
 * at the end of the step.
 * The second one, the ReactionProcess will make the tracks reacting.
 */
class G4ITModelProcessor
{
public:
  G4ITModelProcessor();
  virtual ~G4ITModelProcessor();

  inline void SetModelHandler(G4ITModelHandler*);
  void Initialize();
  //void Initialize(G4MIWorkspace* /*workspace*/);

  void RegisterModel(double time, G4VITStepModel*);

  /**
   * Restore original state of the modelProcessor.
   * This method should be call only by the G4Scheduler
   */
  inline void CleanProcessor();

  G4double CalculateMinTimeStep(G4double currentGlobalTime,
                                G4double definedMinTimeStep);

  void ComputeTrackReaction(G4ITStepStatus fITStepStatus,
                            G4double fGlobalTime,
                            G4double currentTimeStep,
                            G4double previousTimeStep,
                            G4bool reachedUserTimeLimit,
                            G4double fTimeTolerance,
                            G4UserTimeStepAction* fpUserTimeStepAction,
                            G4int fVerbose);

  //____________________________________________________________
  // Time stepper part
  void InitializeStepper(G4double currentGlobalTime,
                         G4double userMinTime);

  bool GetComputeTimeStep()
  {
    return fComputeTimeStep;
  }

protected:

  inline void SetTrack(const G4Track*);

public:
  void CalculateTimeStep(const G4Track*, const G4double);
  void DoCalculateStep();

  //____________________________________________________________
  // Reaction process part
  void FindReaction(G4ITReactionSet* reactionSet,
                    const double currentStepTime,
                    const double previousStepTime,
                    const bool reachedUserStepTimeLimit);

  //____________________________________________________________
  // Get results
  inline const std::vector<std::vector<G4VITStepModel*> >* GetCurrentModel();

  inline std::vector<G4ITReactionChange*>* GetReactionInfo()
  {
    return &fReactionInfo;
  }

  const G4Track* GetTrack() const
  {
    return fpTrack;
  }

  void SetTrackingManager(G4ITTrackingManager* trackingManager)
  {
    fpTrackingManager = trackingManager;
  }

protected:
  void ExtractTimeStepperData();

  /** Copy constructor
   *  \param other Object to copy from
   */
  G4ITModelProcessor(const G4ITModelProcessor& other);
  /** Assignment operator
   *  \param other Object to assign from
   *  \return A reference to this
   */
  G4ITModelProcessor& operator=(const G4ITModelProcessor& other);

  //_____________________________
  // Members

  G4double fTSTimeStep;
  G4ITReactionSet* fReactionSet;
  G4ITTrackingManager* fpTrackingManager;
  G4ITTrackHolder* fpTrackContainer;

  G4bool fInitialized;
  G4ITModelHandler* fpModelHandler;

  const G4Track* fpTrack;
  G4double fUserMinTimeStep;

  // Interactions between different IT types
  // eg : electron/proton
  std::vector<std::vector<G4VITStepModel*> > fCurrentModel;

  // Interactions between ITs of the same type
  // eg : molecule/molecule or electron/electron
  G4VITStepModel* fpModel;
  G4ITModelManager* fpModelManager;

  G4ITType fCurrentType1;
  G4ITType fCurrentType2;

  std::vector<G4ITReactionChange*> fReactionInfo;

  bool fComputeTimeStep;
  bool fComputeReaction;
};

///
// Inline methods
///

inline void G4ITModelProcessor::SetTrack(const G4Track* track)
{
  fpTrack = track;
}

inline const std::vector<std::vector<G4VITStepModel*> >* G4ITModelProcessor::GetCurrentModel()
{
  return &fCurrentModel;
}

inline void G4ITModelProcessor::SetModelHandler(G4ITModelHandler* modelHandler)
{
  if (fInitialized == 1)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "You are trying to set a new model while the model processor has alreaday be initialized";
    G4Exception("G4ITModelProcessor::SetModelHandler", "ITModelProcessor001",
                FatalErrorInArgument, exceptionDescription);
  }
  fpModelHandler = modelHandler;
}

inline void G4ITModelProcessor::CleanProcessor()
{
  fpTrack = 0;
}
#endif // G4ITMODELPROCESSOR_H
