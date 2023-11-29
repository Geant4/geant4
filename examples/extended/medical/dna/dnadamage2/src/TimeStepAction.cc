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
// This example is provided by the Geant4-DNA collaboration
// dnadamage3 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
//
/// \file TimeStepAction.hh
/// \brief Implementation of the TimeStepAction class
///  Filtering of molecules by material at the beggining of the chemical stage

#include "TimeStepAction.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

TimeStepAction::TimeStepAction() : G4UserTimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

TimeStepAction::TimeStepAction(const TimeStepAction& other) :
        G4UserTimeStepAction(other)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

TimeStepAction&
TimeStepAction::operator=(const TimeStepAction& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void TimeStepAction::UserPreTimeStepAction()
{
    G4int Steps = G4Scheduler::Instance()->GetNbSteps();

    if (Steps == 1) {
        G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
        G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
        G4ManyFastLists<G4Track>::iterator it_end   = trackList->end();

        for(;it_begin!=it_end;++it_begin) {
            G4String Material = it_begin->GetStep()->
                         GetPreStepPoint()->GetMaterial()->GetName();
            G4String MoleculeName = GetMolecule(*it_begin)->GetName();
            if (Material != "G4_WATER" && !G4StrUtil::contains(MoleculeName,"DNA")) {
                it_begin->SetTrackStatus(fStopAndKill);
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TimeStepAction::UserPostTimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TimeStepAction::UserReactionAction(const G4Track& /*trackA*/,
    const G4Track& /*trackB*/,
    const std::vector<G4Track*>* /*products*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
