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
// $Id: G4VITProcess.hh 100802 2016-11-02 14:55:27Z gcosmo $
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

#ifndef G4VITProcess_H
#define G4VITProcess_H

#include <G4VProcess.hh>
#include "AddClone_def.hh"
#include "G4ReferenceCast.hh"
#include "G4memory.hh"
#include <typeinfo>

class G4IT;
class G4TrackingInformation;

struct G4ProcessState_Lock
{
  inline virtual ~G4ProcessState_Lock()
  {
    ;
  }
};

/*
 class G4ProcessStateHandle_Lock : public G4shared_ptr<G4ProcessState_Lock>
 {
 public:
 G4ProcessStateHandle_Lock(G4ProcessState_Lock* plock) : G4shared_ptr<G4ProcessState_Lock>(plock)
 {}
 virtual ~G4ProcessStateHandle_Lock(){}
 };
 */

#define InitProcessState(destinationType,source) \
    reference_cast<destinationType>(source)

#define DowncastProcessState(destinationType) \
	G4dynamic_pointer_cast<destinationType>(G4VITProcess::fpState)

#define UpcastProcessState(destinationType) \
	G4dynamic_pointer_cast<destinationType>(G4VITProcess::fpState)

#define DowncastState(destinationType,source) \
	G4dynamic_pointer_cast<destinationType>(source)

#define UpcastState(destinationType,source) \
	G4dynamic_pointer_cast<destinationType>(source)

/**
 * G4VITProcess inherits from G4VProcess.
 * A G4VITProcess is able to save its current state for a given track into G4IT.
 * This state may be retrieve latter on to be used by the G4VITProcess.
 * Each G4VITProcess is tagged.
 */

class G4VITProcess : public G4VProcess
{
public:
  //__________________________________
  // Constructors & destructors
  G4VITProcess(const G4String& name, G4ProcessType type = fNotDefined);

  virtual ~G4VITProcess();
  G4VITProcess(const G4VITProcess& other);
  G4VITProcess& operator=(const G4VITProcess& other);

  // equal opperators
  G4int operator==(const G4VITProcess &right) const;
  G4int operator!=(const G4VITProcess &right) const;

  G4IT_TO_BE_CLONED(G4VITProcess)

  size_t GetProcessID() const
  {
    return fProcessID;
  }

//    G4ProcessState_Lock* GetProcessState()
//    {
//        return fpState;
//    }
//
//    void SetProcessState(G4ProcessState_Lock* aProcInfo)
//    {
//        fpState = (G4ProcessState*) aProcInfo;
//    }

  G4shared_ptr<G4ProcessState_Lock> GetProcessState()
  {
    return UpcastProcessState(G4ProcessState_Lock);
  }

  void SetProcessState(G4shared_ptr<G4ProcessState_Lock> aProcInfo)
  {
    fpState = DowncastState(G4ProcessState, aProcInfo);
  }

  void ResetProcessState()
  {
    fpState.reset();
  }

  //__________________________________
  // Initialize and Save process info

  virtual void StartTracking(G4Track*);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&)
  {
  }

  inline G4double GetInteractionTimeLeft();

  /** WARNING : Redefine the method of G4VProcess
   * reset (determine the value of)NumberOfInteractionLengthLeft
   */
  virtual void ResetNumberOfInteractionLengthLeft();

  inline G4bool ProposesTimeStep() const;

  inline static const size_t& GetMaxProcessIndex();

protected:
  // with description

  void RetrieveProcessInfo();
  void CreateInfo();

  //__________________________________
  // Process info
  // friend class G4TrackingInformation ;

  struct G4ProcessState : public G4ProcessState_Lock
  {
  public:
    G4ProcessState();
    virtual ~G4ProcessState();

    virtual G4String GetType()
    {
      return "G4ProcessState";
    }

    G4double theNumberOfInteractionLengthLeft;
    // The flight length left for the current tracking particle
    // in unit of "Interaction length".

    G4double theInteractionTimeLeft;
    // Time left before the interaction : for at rest processes

    G4double currentInteractionLength;
    // The InteractionLength in the current material

    template<typename T>
      T* GetState()
      {
        return dynamic_cast<T*>(this);
      }
  };

  template<typename T>
    class G4ProcessStateBase : public G4ProcessState
    {
    public:
      G4ProcessStateBase() :
          G4ProcessState()
      {
      }
      virtual ~G4ProcessStateBase()
      {
      }

      virtual G4String GetType()
      {
        return typeid(T).name();
      }
    };

  template<typename T>
    T* GetState()
    {
      return fpState->GetState<T>();
    }

  G4shared_ptr<G4ProcessState> fpState;

  void virtual SubtractNumberOfInteractionLengthLeft(G4double previousStepSize);

  inline virtual void ClearInteractionTimeLeft();

  inline virtual void ClearNumberOfInteractionLengthLeft();
  // clear NumberOfInteractionLengthLeft
  // !!! This method should be at the end of PostStepDoIt()
  // !!! and AtRestDoIt
  //_________________________________________________

  void SetInstantiateProcessState(G4bool flag)
  {
    fInstantiateProcessState = flag;
  }

  G4bool InstantiateProcessState()
  {
    return fInstantiateProcessState;
  }

  G4bool fProposesTimeStep;

private:

  size_t fProcessID;
  // During all the simulation will identify a process, so if two identical
  // processes are created using a copy constructor they will have the same
  // fProcessID. NOTE: due to MT, this cannot be "const".

  static/*G4ThreadLocal*/size_t *fNbProcess;

  G4bool fInstantiateProcessState;
  //_________________________________________________
  // Redefine needed members and method of G4VProcess
  G4double* theNumberOfInteractionLengthLeft;
  G4double* currentInteractionLength;
  G4double* theInteractionTimeLeft;
};

inline void G4VITProcess::ClearInteractionTimeLeft()
{
  fpState->theInteractionTimeLeft = -1.0;
}

inline void G4VITProcess::ClearNumberOfInteractionLengthLeft()
{
  fpState->theNumberOfInteractionLengthLeft = -1.0;
}

inline void G4VITProcess::ResetNumberOfInteractionLengthLeft()
{
  fpState->theNumberOfInteractionLengthLeft = -std::log( G4UniformRand());
}

inline G4double G4VITProcess::GetInteractionTimeLeft()
{
  if (fpState) return fpState->theInteractionTimeLeft;

  return -1;
}

inline G4bool G4VITProcess::ProposesTimeStep() const
{
  return fProposesTimeStep;
}

inline const size_t& G4VITProcess::GetMaxProcessIndex()
{
  if (!fNbProcess) fNbProcess = new size_t(0);
  return *fNbProcess;
}

inline
void G4VITProcess::SubtractNumberOfInteractionLengthLeft(G4double previousStepSize)
{
  if (fpState->currentInteractionLength > 0.0)
  {
    fpState->theNumberOfInteractionLengthLeft -= previousStepSize
        / fpState->currentInteractionLength;
    if (fpState->theNumberOfInteractionLengthLeft < 0.)
    {
      fpState->theNumberOfInteractionLengthLeft = CLHEP::perMillion;
    }

  }
  else
  {
#ifdef G4VERBOSE
    if (verboseLevel > 0)
    {
      G4cerr << "G4VITProcess::SubtractNumberOfInteractionLengthLeft()";
      G4cerr << " [" << theProcessName << "]" << G4endl;
      G4cerr << " currentInteractionLength = "
             << fpState->currentInteractionLength << " [mm]";
      G4cerr << " previousStepSize = " << previousStepSize << " [mm]";
      G4cerr << G4endl;
    }
#endif
             G4String msg = "Negative currentInteractionLength for ";
             msg += theProcessName;
             G4Exception("G4VITProcess::SubtractNumberOfInteractionLengthLeft()",
                "ProcMan201",EventMustBeAborted,
                msg);
  }
}
#endif // G4VITProcess_H
