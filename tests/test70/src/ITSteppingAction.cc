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
#include "ITSteppingAction.hh"
// #include "TTree.h"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"
#include "G4UnitsTable.hh"

ITSteppingAction::ITSteppingAction() : G4UserSteppingAction()
{

}

ITSteppingAction::~ITSteppingAction()
{
}

void ITSteppingAction::UserSteppingAction(const G4Step* )
{
    /*
     *    G4Track* track = step->GetTrack();
     *    G4cout << "Track ID : " << track->GetTrackID()
     *      << "\t"
     *      << "Molecule name : " << GetMolecule(track)->GetName()
     *      << "\t"
     *      <<"Track time " << G4BestUnit(track->GetGlobalTime(), "Time")
     *      << " stepping " << G4endl;
     */

    /*
    //    G4cout << "ITSteppingAction::UserSteppingAction" << G4endl;
    G4Track* track = step->GetTrack();
    G4Molecule* molecule = GetMolecule(track);

    TTree* tree = 0;

    const G4String& moleculeName = molecule->GetName();

    fMolName = moleculeName.c_str();
    fTrackID= track->GetTrackID() ;
    fTime = track->GetGlobalTime() /picosecond;
    const G4ThreeVector& positionUm = track->GetPosition() /micrometer;
    fPosition.SetXYZ(positionUm.getX(),positionUm.getY(),positionUm.getZ());

    map<const G4Molecule, TTree*>::iterator it = fpTree.find(*molecule);
    if(it == fpTree.end())
    {
        if(track->GetParticleDefinition() == G4H2O::Definition()) return;
        // put there in order to avoid checking at every step
        // you probably don't want to register the water molecule position
        // if you do, you should also record their electronic configuration
        // (eg.= excited : level, ionised : level)

        G4String molPos = "molPos_" + moleculeName;
        tree = new TTree(molPos,molPos);
        fpTree[*molecule]  = tree;
        tree->Branch("time",&fTime,"time/D");
        tree->Branch("trackID",&fTrackID,"trackID/I");
        tree->Branch("position",&fPosition);
        tree->Branch("name",(void*)fMolName,"name/C");
    }
    else
    {
        tree = it->second;
    }

    tree->Fill();

    if(track->GetCurrentStepNumber() == 1) // To save the first position/time
    {
        G4StepPoint* preStepPoint = step->GetPreStepPoint();
        fTime = preStepPoint->GetGlobalTime() /picosecond;
        const G4ThreeVector& positionUm2 = preStepPoint->GetPosition() /micrometer;
        fPosition.SetXYZ(positionUm2.getX(),positionUm2.getY(),positionUm2.getZ());
        tree->Fill();
    }
    */
}
