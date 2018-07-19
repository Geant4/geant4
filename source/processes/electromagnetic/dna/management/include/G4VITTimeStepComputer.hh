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
// $Id: G4VITTimeStepComputer.hh 100802 2016-11-02 14:55:27Z gcosmo $
////
// Author: Mathieu Karamitros
////
// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, so do not hesitate to
// send us your feedback!
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial.
// If you use Geant4-DNA chemistry and you publish papers about
// your software, in addition to the general paper on Geant4-DNA:
//
// The Geant4-DNA project, S. Incerti et al.,
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following papers reference papers
// related to chemistry:
//
// Diffusion-controlled reactions modelling in Geant4-DNA,
// M. Karamitros et al., 2014 (submitted)
// Modeling Radiation Chemistry in the Geant4 Toolkit, M. Karamitros et al.,
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508

#ifndef G4VITTimeStepper_H
#define G4VITTimeStepper_H

#include "G4Track.hh"
#include "G4ITReactionTable.hh"
#include "G4ReferenceCountedHandle.hh"
#include "AddClone_def.hh"
#include "G4memory.hh"

//typedef G4ReferenceCountedHandle< std::vector<G4Track*> > G4TrackVectorHandle;
typedef G4shared_ptr< std::vector<G4Track*> > G4TrackVectorHandle;

/**
  * Before stepping all tracks G4Scheduler calls all the G4VITModel
  * which may contain a G4VITTimeStepper (optionnal).
  * G4VITTimeStepper returns what should be the next global time step.
  * Time step that will be used to step all tracks.
  */

class G4VITTimeStepComputer
{
public:
    G4VITTimeStepComputer();
    virtual ~G4VITTimeStepComputer();

    G4VITTimeStepComputer(const G4VITTimeStepComputer&);
    G4VITTimeStepComputer& operator=(const G4VITTimeStepComputer& other);

    /** This macro defined in AddClone_def **/
    G4IT_TO_BE_CLONED(G4VITTimeStepComputer)

    // First initialization (done once for all at the begin of the run)
    // eg. check if the reaction table is given ...
    inline virtual void Initialize(){;}

    // Preparation part
    static void SetTimes(const G4double&, const G4double&);
//    inline virtual void PrepareForAllProcessors(){;}
    inline virtual void Prepare() ;

    virtual G4double CalculateStep(const G4Track&, const G4double&) = 0;

    inline G4TrackVectorHandle GetReactants();
    inline virtual void ResetReactants()
//    {fReactants = 0;}
    {fReactants.reset();}

    //
    inline G4double GetSampledMinTimeStep() ;
    
    inline void SetReactionTable(const G4ITReactionTable*);
    inline const G4ITReactionTable* GetReactionTable();

protected :
    static G4ThreadLocal G4double fCurrentGlobalTime ;
    static G4ThreadLocal G4double fUserMinTimeStep   ;

    G4double fSampledMinTimeStep ;
    G4TrackVectorHandle fReactants;

    const G4ITReactionTable* fpReactionTable;

private:
    G4int fVerbose ;
};

inline void G4VITTimeStepComputer::SetReactionTable(const G4ITReactionTable* table)
{
    fpReactionTable = table;
}

inline const G4ITReactionTable* G4VITTimeStepComputer::GetReactionTable()
{
    return fpReactionTable ;
}

inline void G4VITTimeStepComputer::Prepare()
{
//    fReactants = 0 ;
	fReactants.reset() ;
}

inline G4double G4VITTimeStepComputer::GetSampledMinTimeStep()
{
    return fSampledMinTimeStep ;
}

inline G4TrackVectorHandle G4VITTimeStepComputer::GetReactants()
{
    return  fReactants ;
}
#endif // G4VITTimeStepper_H
