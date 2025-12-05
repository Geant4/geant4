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
/// \file TimeStepAction.hh
/// \brief Definition of the TimeStepAction class

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//

#ifndef ITACTION_H
#define ITACTION_H

#include "G4UserTimeStepAction.hh"
class G4MolecularConfiguration;

class TimeStepAction : public G4UserTimeStepAction
{
  public:
    TimeStepAction();
    ~TimeStepAction() override = default;
    TimeStepAction(const TimeStepAction& other) = delete;
    TimeStepAction& operator=(const TimeStepAction& other) = delete;

    void StartProcessing() override { ; }

    static void PrintSpecieInfo(G4MolecularConfiguration* molconf);

    /** In this method, the user can use :
     * G4ITTimeStepper::Instance()->GetGlobalTime(),
     *    to know the current simulation time
     * G4ITTimeStepper::Instance()->GetTimeStep(),
     *    to know the selected minimum time
     * WARNING :
     *    The call of this method happens before the call of DoIT methods
     */
    void UserPreTimeStepAction() override { ; }
    void UserPostTimeStepAction() override;

    /**
     * Inform about a reaction
     */
    void UserReactionAction(const G4Track& /*trackA*/, const G4Track& /*trackB*/,
                                    const std::vector<G4Track*>* /*products*/) override;

    void EndProcessing() override { ; }
};

#endif  // ITACTION_H
