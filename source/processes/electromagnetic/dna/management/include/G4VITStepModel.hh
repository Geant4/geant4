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
// $Id: G4VITStepModel.hh 100802 2016-11-02 14:55:27Z gcosmo $
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

#ifndef G4VITMODEL_HH
#define G4VITMODEL_HH

#include "AddClone_def.hh"
#include "G4VITTimeStepComputer.hh"
#include "G4VITReactionProcess.hh"
#include "G4ITReactionTable.hh"

/**
 * Define what to do before stepping and after stepping.
 * The concrete implementation of G4VITModel defines the interaction
 * between two G4IT types. The type might be just equal like :
 * Molecule + Molecule, or different : Molecule + Atom.
 */

class G4VITStepModel
{
public:
    /** Default constructor */
    G4VITStepModel(const G4String& aName = "NoName");
    /** Default destructor */
    virtual ~G4VITStepModel();

    /* Macro define in AddClone_def*/
    G4IT_TO_BE_CLONED(G4VITStepModel)

    void IsApplicable(G4ITType& type1, G4ITType& type2) ;
    void virtual PrintInfo(){;}

    virtual void Initialize();

    inline void SetTimeStepper(G4VITTimeStepComputer* timeStepper);
    inline void SetReactionProcess(G4VITReactionProcess* reactionProcess);

    inline G4VITTimeStepComputer* GetTimeStepper();
    inline const G4String& GetName();

    inline G4VITReactionProcess* GetReactionProcess();
    inline void SetReactionTable(G4ITReactionTable*);
    inline const G4ITReactionTable* GetReactionTable();

protected:

    G4String fName;

    G4VITTimeStepComputer* fpTimeStepper;
    G4VITReactionProcess* fpReactionProcess;

    const G4ITReactionTable* fpReactionTable ;

    G4ITType fType1;
    G4ITType fType2;

protected :
    /** Copy constructor
         *  \param other Object to copy from
         */
    G4VITStepModel(const G4VITStepModel& other);
    /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
    G4VITStepModel& operator=(const G4VITStepModel& other);
};

inline void G4VITStepModel::SetReactionTable(G4ITReactionTable* table)
{
    fpReactionTable = table;
}

inline const G4ITReactionTable* G4VITStepModel::GetReactionTable()
{
    return fpReactionTable ;
}

inline void G4VITStepModel::SetTimeStepper(G4VITTimeStepComputer* timeStepper)
{
    fpTimeStepper = timeStepper ;
}

inline void G4VITStepModel::SetReactionProcess(G4VITReactionProcess* reactionProcess)
{
    fpReactionProcess = reactionProcess ;
}

inline G4VITTimeStepComputer* G4VITStepModel::GetTimeStepper()
{
    return fpTimeStepper;
}

inline G4VITReactionProcess* G4VITStepModel::GetReactionProcess()
{
    return fpReactionProcess ;
}

inline const G4String& G4VITStepModel::GetName()
{
    return fName;
}

#endif // G4VITMODEL_HH
