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
// $Id: G4VITReactionProcess.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4VITREACTIONPROCESS_H
#define G4VITREACTIONPROCESS_H

#include "globals.hh"
#include "G4IT.hh"
#include "G4ITReactionChange.hh"
#include "G4ITReactionTable.hh"
#include "AddClone_def.hh"

/**
 * G4VITReactionProcess defines the reaction between two G4IT.
 * It should be stored in a G4VITModel.
 */

class G4VITReactionProcess
{
public:
    /** Default constructor */
    G4VITReactionProcess();
    /** Default destructor */
    virtual ~G4VITReactionProcess();
    /** Copy constructor
         *  \param other Object to copy from
         */
    G4VITReactionProcess(const G4VITReactionProcess& other);
    /** Will Clone the reaction process
         * i.e. new reaction process will be created with same features
         * as the parent one.
         * Use preprocessor macro : AddCloneReactionProcess
         * the copy constructor of the derived class must be
         * implemented.
         * This macro is defined in AddClone_def
         */
    G4IT_TO_BE_CLONED(G4VITReactionProcess)

    /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
    G4VITReactionProcess& operator=(const G4VITReactionProcess& other);

    /** First initialization (done once for all at the begin of the run)
     * eg. check if the reaction table is given ...
     */
    virtual void Initialize(){;}

    virtual G4bool IsApplicable(G4ITType,G4ITType) const {return true;}
//    virtual void GetApplicableTypes(G4ITType&, G4ITType&) const;

    virtual G4bool TestReactibility(const G4Track&,
                                    const G4Track&,
                                    const double /*currentStepTime*/,
                                    const double /*previousStepTime*/,
                                    bool /*reachedUserStepTimeLimit*/) = 0;

    virtual G4ITReactionChange* MakeReaction(const G4Track&, const G4Track&) = 0;

    inline void SetReactionTable(const G4ITReactionTable*);

    inline void ResetChanges();

protected:
//    G4ITType    fApplicableType1;
//    G4ITType    fApplicableType2;

    const G4ITReactionTable* fpReactionTable;
    G4ITReactionChange*  fpChanges ;
    G4String    fName ;
};

inline void G4VITReactionProcess::SetReactionTable(const G4ITReactionTable* table)
{
    fpReactionTable = table;
}

inline void G4VITReactionProcess::ResetChanges()
{
    fpChanges = 0;
}
#endif // G4VITREACTIONPROCESS_H
