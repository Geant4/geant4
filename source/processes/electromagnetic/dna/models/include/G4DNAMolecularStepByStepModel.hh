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
// $Id: G4DNAMolecularStepByStepModel.hh 94218 2015-11-09 08:24:48Z gcosmo $
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


#ifndef G4MOLECULARSTEPBYSTEPMODEL_H
#define G4MOLECULARSTEPBYSTEPMODEL_H

#include <AddClone_def.hh>
#include <G4String.hh>
#include <G4VITStepModel.hh>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;


/**
  * G4DNAMolecularStepByStepModel :
  *  - TimeStepper : G4DNAMolecularEncounterStepper
  *  - ReactionProcess : G4DNAMolecularReaction
  * Before each step, the next minimum encounter time is calculated for each
  * pair of molecule. The minimum time step is selected. All the molecules are stepped
  * within this time step. Then, only the relevant pair of molecules are checked for
  * reaction.
  */

class G4DNAMolecularStepByStepModel : public G4VITStepModel
{
public:
    /** Default constructor */
    G4DNAMolecularStepByStepModel(const G4String& name = "DNAMolecularStepByStepModel");
    /** Default destructor */
    virtual ~G4DNAMolecularStepByStepModel();

    G4DNAMolecularStepByStepModel(const G4DNAMolecularStepByStepModel&);

    G4IT_ADD_CLONE(G4VITStepModel, G4DNAMolecularStepByStepModel)

    virtual void PrintInfo();
    virtual void Initialize();

    inline void SetReactionModel(G4VDNAReactionModel*);
    inline G4VDNAReactionModel* GetReactionModel();

protected:
    const G4DNAMolecularReactionTable*& fMolecularReactionTable;
    G4VDNAReactionModel* fReactionModel;

private :
    /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
    G4DNAMolecularStepByStepModel& operator=(const G4DNAMolecularStepByStepModel& /*other*/);
};

inline void G4DNAMolecularStepByStepModel::SetReactionModel(G4VDNAReactionModel* reactionModel)
{
    fReactionModel = reactionModel;
}

inline G4VDNAReactionModel* G4DNAMolecularStepByStepModel::GetReactionModel()
{
    return fReactionModel;
}

#endif // G4MOLECULARSTEPBYSTEPMODEL_H
