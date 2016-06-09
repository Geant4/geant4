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
// $Id: G4VITTimeStepper.hh 64057 2012-10-30 15:04:49Z gcosmo $
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
#ifndef G4VITTimeStepper_H
#define G4VITTimeStepper_H

#include "G4Track.hh"
#include "G4ITReactionTable.hh"
#include "G4ReferenceCountedHandle.hh"
#include "AddClone_def.hh"

typedef G4ReferenceCountedHandle< std::vector<G4Track*> > G4TrackVectorHandle;

/**
  * Before stepping all tracks G4ITStepManager calls all the G4VITModel
  * which may contain a G4VITTimeStepper (optionnal).
  * G4VITTimeStepper returns what should be the next global time step.
  * Time step that will be used to step all tracks.
  */

class G4VITTimeStepper
{
public:
    G4VITTimeStepper();
    virtual ~G4VITTimeStepper();

    G4VITTimeStepper(const G4VITTimeStepper&);
    G4VITTimeStepper& operator=(const G4VITTimeStepper& other);

    /** This macro defined in AddClone_def **/
    G4IT_TO_BE_CLONED(G4VITTimeStepper)

    // First initialization (done once for all at the begin of the run)
    // eg. check if the reaction table is given ...
    inline virtual void Initialize(){;}

    // Preparation part
    static void SetTimes(const G4double&, const G4double&);
//    inline virtual void PrepareForAllProcessors(){;}
    inline virtual void Prepare() ;

    virtual G4double CalculateStep(const G4Track&, const G4double&) = 0;

    inline G4TrackVectorHandle GetReactants();
    inline virtual void ResetReactants(){fReactants = 0;}

    //
    inline G4double GetSampledMinTimeStep() ;
    
    inline void SetReactionTable(const G4ITReactionTable*);
    inline const G4ITReactionTable* GetReactionTable();

protected :
    static G4double fCurrentGlobalTime ;
    static G4double fUserMinTimeStep   ;

    G4double fSampledMinTimeStep ;
    G4TrackVectorHandle fReactants;

    const G4ITReactionTable* fpReactionTable;

private:
    G4int fVerbose ;
};

inline void G4VITTimeStepper::SetReactionTable(const G4ITReactionTable* table)
{
    fpReactionTable = table;
}

inline const G4ITReactionTable* G4VITTimeStepper::GetReactionTable()
{
    return fpReactionTable ;
}

inline void G4VITTimeStepper::Prepare()
{
    fReactants = 0 ;
}

inline G4double G4VITTimeStepper::GetSampledMinTimeStep()
{
    return fSampledMinTimeStep ;
}

inline G4TrackVectorHandle G4VITTimeStepper::GetReactants()
{
    return  fReactants ;
}
#endif // G4VITTimeStepper_H
