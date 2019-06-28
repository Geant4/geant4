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


#pragma once

#include "globals.hh"
#include <memory>

class G4ITReactionTable;
class G4ITReactionChange;
class G4Track;
struct G4ITType;

/**
 * G4VITReactionProcess defines the reaction between two G4IT.
 * It should be stored in a G4VITModel.
 * \deprecated This class will be removed
 */
class G4VITReactionProcess
{
public:
    G4VITReactionProcess() = default;
    virtual ~G4VITReactionProcess() = default;
    G4VITReactionProcess(const G4VITReactionProcess& other) = delete;
    G4VITReactionProcess& operator=(const G4VITReactionProcess& other) = delete;

    /** First initialization (done once for all at the begin of the run)
     * eg. check if the reaction table is given ...
     */
    virtual void Initialize(){;}

    virtual G4bool IsApplicable(const G4ITType&, const G4ITType&) const;

    virtual G4bool TestReactibility(const G4Track&,
                                    const G4Track&,
                                    double /*currentStepTime*/,
                                    bool /*reachedUserStepTimeLimit*/) = 0;

    virtual std::unique_ptr<G4ITReactionChange> MakeReaction(const G4Track&, const G4Track&) = 0;

    virtual void SetReactionTable(const G4ITReactionTable*);

protected:
    const G4ITReactionTable* fpReactionTable = nullptr;
    G4String fName ;
};
