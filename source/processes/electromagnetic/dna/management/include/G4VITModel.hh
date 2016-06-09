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
#ifndef G4VITMODEL_HH
#define G4VITMODEL_HH

#include "AddClone_def.hh"
#include "G4VITTimeStepper.hh"
#include "G4VITReactionProcess.hh"
#include "G4ITReactionTable.hh"

/**
 * Define what to do before stepping and after stepping.
 * The concrete implementation of G4VITModel defines the interaction
 * between two G4IT types. The type might be just equal like :
 * Molecule + Molecule, or different : Molecule + Atom.
 */

class G4VITModel
{
public:
    /** Default constructor */
    G4VITModel(const G4String& aName = "NoName");
    /** Default destructor */
    virtual ~G4VITModel();

    /* Macro define in AddClone_def*/
    G4IT_TO_BE_CLONED(G4VITModel)

    void IsApplicable(G4ITType& type1, G4ITType& type2) ;
    void virtual PrintInfo(){;}

    virtual void Initialize();

    inline void SetTimeStepper(G4VITTimeStepper* timeStepper);
    inline void SetReactionProcess(G4VITReactionProcess* reactionProcess);

    inline G4VITTimeStepper* GetTimeStepper();
    inline const G4String& GetName();

    inline G4VITReactionProcess* GetReactionProcess();
    inline void SetReactionTable(G4ITReactionTable*);
    inline const G4ITReactionTable* GetReactionTable();

protected:

    G4String fName;

    G4VITTimeStepper* fTimeStepper;
    G4VITReactionProcess* fReactionProcess;

    const G4ITReactionTable* fReactionTable ;

    G4ITType fType1;
    G4ITType fType2;

protected :
    /** Copy constructor
         *  \param other Object to copy from
         */
    G4VITModel(const G4VITModel& other);
    /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
    G4VITModel& operator=(const G4VITModel& other);
};

inline void G4VITModel::SetReactionTable(G4ITReactionTable* table)
{
    fReactionTable = table;
}

inline const G4ITReactionTable* G4VITModel::GetReactionTable()
{
    return fReactionTable ;
}

inline void G4VITModel::SetTimeStepper(G4VITTimeStepper* timeStepper)
{
    fTimeStepper = timeStepper ;
}

inline void G4VITModel::SetReactionProcess(G4VITReactionProcess* reactionProcess)
{
    fReactionProcess = reactionProcess ;
}

inline G4VITTimeStepper* G4VITModel::GetTimeStepper()
{
    return fTimeStepper;
}

inline G4VITReactionProcess* G4VITModel::GetReactionProcess()
{
    return fReactionProcess ;
}

inline const G4String& G4VITModel::GetName()
{
    return fName;
}

#endif // G4VITMODEL_HH
