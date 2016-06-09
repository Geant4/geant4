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
// $Id: G4VITProcess.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4VITProcess_H
#define G4VITProcess_H

#include <G4VProcess.hh>
#include "AddClone_def.hh"
#include "G4ReferenceCast.hh"

class G4IT ;
class G4TrackingInformation ;

struct G4ProcessState_Lock{
    inline virtual ~G4ProcessState_Lock(){;}
};

#define InitProcessState(destination,source) \
    destination(reference_cast(destination,source))

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

    G4ProcessState_Lock* GetProcessState()
    {
        return fpState;
    }

    void SetProcessState(G4ProcessState_Lock* aProcInfo)
    {
        fpState = (G4ProcessState*) aProcInfo;
    }

    //__________________________________
    // Initialize and Save process info

    virtual void StartTracking(G4Track*);

    virtual void BuildPhysicsTable(const G4ParticleDefinition&){}

    inline G4double GetInteractionTimeLeft();

    /** WARNING : Redefine the method of G4VProcess
    * reset (determine the value of)NumberOfInteractionLengthLeft
    */
    virtual void  ResetNumberOfInteractionLengthLeft();

    inline G4bool ProposesTimeStep() const;

    inline static const size_t& GetMaxProcessIndex();

protected:  // with description

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

        G4double          theNumberOfInteractionLengthLeft;
        // The flight length left for the current tracking particle
        // in unit of "Interaction length".

        G4double          theInteractionTimeLeft;
        // Time left before the interaction : for at rest processes

        G4double          currentInteractionLength;
        // The InteractionLength in the current material
    };

    G4ProcessState* fpState ;

    inline virtual void ClearInteractionTimeLeft();

    //_________________________________________________
    // Redefine needed members and method of G4VProcess
    virtual void      SubtractNumberOfInteractionLengthLeft(
        G4double previousStepSize
    );
    // subtract NumberOfInteractionLengthLeft by the value corresponding to
    // previousStepSize

    inline virtual void      ClearNumberOfInteractionLengthLeft();
    // clear NumberOfInteractionLengthLeft
    // !!! This method should be at the end of PostStepDoIt()
    // !!! and AtRestDoIt
    //_________________________________________________


    void SetInstantiateProcessState(G4bool flag)
    { fInstantiateProcessState = flag; }

    G4bool InstantiateProcessState() { return fInstantiateProcessState; }

    G4bool fProposesTimeStep;

private :
    const size_t fProcessID; // During all the simulation will identify a
    // process, so if two identical process are created using a copy constructor
    // they will have the same fProcessID
    static size_t fNbProcess ;

    G4bool fInstantiateProcessState;
    //_________________________________________________
    // Redefine needed members and method of G4VProcess
    G4double*          theNumberOfInteractionLengthLeft;
    G4double*          currentInteractionLength;
    G4double*          theInteractionTimeLeft;
};

inline void G4VITProcess::ClearInteractionTimeLeft()
{
    fpState->theInteractionTimeLeft = -1.0;
}

inline void G4VITProcess::ClearNumberOfInteractionLengthLeft()
{
    fpState->theNumberOfInteractionLengthLeft =  -1.0;
}

inline void G4VITProcess::ResetNumberOfInteractionLengthLeft()
{
    fpState->theNumberOfInteractionLengthLeft =  -std::log( G4UniformRand() );
}

inline G4double G4VITProcess::GetInteractionTimeLeft()
{
    if(fpState)
        return fpState->theInteractionTimeLeft ;

    return -1 ;
}

inline G4bool G4VITProcess::ProposesTimeStep() const
{
    return fProposesTimeStep;
}

inline const size_t& G4VITProcess::GetMaxProcessIndex()
{
    return fNbProcess ;
}
#endif // G4VITProcess_H
