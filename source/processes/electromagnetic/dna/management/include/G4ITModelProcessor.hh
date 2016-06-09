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
// $Id: G4ITModelProcessor.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4ITMODELPROCESSOR_H
#define G4ITMODELPROCESSOR_H

#include <vector>
#include "G4ITReactionChange.hh"
#include "G4ITType.hh"
#include "G4ITModelHandler.hh"

class G4VITTimeStepper;
class G4VITReactionProcess;
class G4ITModelHandler;

typedef G4ReferenceCountedHandle< std::vector<G4Track*> > G4TrackVectorHandle;

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
    /** Default constructor */
    G4ITModelProcessor();
    /** Default destructor */
    virtual ~G4ITModelProcessor();


    inline void SetModelHandler(G4ITModelHandler*);

    void Initialize();
    
    /**
     * Restaure original state of the modelProcessor.
     * This method should be call only by the ITStepManager
     */ 
    inline void CleanProcessor();

    //____________________________________________________________
    // Time stepper part
    void InitializeStepper(const G4double& currentGlobalTime,
                           const G4double& userMinTime);
protected :

    inline void SetTrack(const G4Track*);

public :
    void CalculateTimeStep(const G4Track*, const G4double);
    void DoCalculateStep();

    //____________________________________________________________
    // Reaction process part
    void FindReaction(std::map<G4Track*, G4TrackVectorHandle>*,
                      const double currentStepTime,
                      const double previousStepTime,
                      const bool reachedUserStepTimeLimit) ;

    //____________________________________________________________
    // Get results
    inline const std::vector<std::vector<G4VITModel*> >* GetCurrentModel();

    inline std::vector<G4ITReactionChange*>*   GetReactionInfo()
    {
        return &fReactionInfo;
    }

    const G4Track* GetTrack() const
    {
        return fpTrack;
    }


protected:
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
    G4bool fInitialized;
    G4ITModelHandler* fpModelHandler ;

    const G4Track* fpTrack;
    G4double fUserMinTimeStep;

    // Attributes for interaction between many IT types
    // eg : electron/proton
    std::vector<std::vector<G4VITModel*> >   fCurrentModel;

    // Attributes for interaction between one type of IT
    // eg : molecule/molecule or electron/electron
    G4VITModel*                 fpModel;
    G4ITModelManager*           fpModelManager;
//    G4double                    fNextTimeChangeModel ;

    G4ITType    fCurrentType1;
    G4ITType    fCurrentType2;

    // Atribute for reactions
    std::vector<G4ITReactionChange*> fReactionInfo ;
    static std::map<const G4Track*, G4bool> fHasReacted;
};

///
// Inline methods
///

inline void G4ITModelProcessor::SetTrack(const G4Track* track)
{
    fpTrack = track;
}

inline const std::vector<std::vector<G4VITModel*> >* G4ITModelProcessor::GetCurrentModel()
{
    return &fCurrentModel ;
}

inline void G4ITModelProcessor::SetModelHandler(G4ITModelHandler* modelHandler)
{
    if(fInitialized == 1)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "You are trying to set a new model while the model processor has alreaday be initialized";
        G4Exception("G4ITModelProcessor::SetModelHandler","ITModelProcessor001",
                    FatalErrorInArgument,exceptionDescription);
    }
    fpModelHandler = modelHandler;
}

inline void G4ITModelProcessor::CleanProcessor()
{
    fpTrack = 0;
}
#endif // G4ITMODELPROCESSOR_H
