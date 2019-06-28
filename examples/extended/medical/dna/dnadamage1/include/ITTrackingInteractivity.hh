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
/// \file ITTrackingInteractivity.hh
/// \brief Definition of the ITTrackingInteractivity class

#pragma once
#include "G4ITTrackingInteractivity.hh"
#include <vector>
#include <memory>

class G4VTrajectory;

class ITTrackingInteractivity: public G4ITTrackingInteractivity
{
public:
    ITTrackingInteractivity();
    ~ITTrackingInteractivity() override;

    void Initialize() override;
    void StartTracking(G4Track*) override;
    void AppendStep(G4Track* track, G4Step* step) override;
    void EndTracking(G4Track*) override;
    void Finalize() override;

    void SetUserAction(G4UserTrackingAction*);
    inline std::unique_ptr<G4UserTrackingAction> 
    GetUserTrackingAction();

    void SetUserAction(G4UserSteppingAction*);
    inline std::unique_ptr<G4UserSteppingAction> 
    GetUserSteppingAction();

private:
    int fStoreTrajectory;
    std::vector<G4VTrajectory*> fTrajectories;
    std::unique_ptr<G4UserTrackingAction> 
    fpUserTrackingAction;
    std::unique_ptr<G4UserSteppingAction> 
    fpUserSteppingAction;
};

inline
void ITTrackingInteractivity::SetUserAction(
    G4UserTrackingAction* trackAct)
{
    fpUserTrackingAction.reset(trackAct);
}

inline
void ITTrackingInteractivity::SetUserAction(
    G4UserSteppingAction* stepAct)
{
    fpUserSteppingAction.reset(stepAct);
}

inline std::unique_ptr<G4UserSteppingAction>
ITTrackingInteractivity::GetUserSteppingAction()
{
    return std::move(fpUserSteppingAction);
}

inline std::unique_ptr<G4UserTrackingAction>
ITTrackingInteractivity::GetUserTrackingAction()
{
    return std::move(fpUserTrackingAction);
}

