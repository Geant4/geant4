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
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file TimeStepAction.hh
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4IT.hh"

#include "G4ITTrackHolder.hh"
#include "G4Molecule.hh"
#include "G4Scheduler.hh"
#include "G4AnalysisManager.hh"

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
//using namespace G4DNAPARSER;

TimeStepAction::TimeStepAction() :
    G4UserTimeStepAction(),
    fpDetector(0)
{
  fpDetector =
      dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
          ->GetUserDetectorConstruction());

  AddTimeStep(1     * picosecond, 0.1 * picosecond);
  AddTimeStep(10    * picosecond, 1   * picosecond);
  AddTimeStep(100   * picosecond, 3   * picosecond);
  AddTimeStep(1000  * picosecond, 10  * picosecond);
  AddTimeStep(10000 * picosecond, 100 * picosecond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::~TimeStepAction()
{
  //dtor
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction(const TimeStepAction& other) :
    G4UserTimeStepAction(other)
{
  //copy ctor
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction&
TimeStepAction::operator=(const TimeStepAction& rhs)
{
  if (this == &rhs)
    return *this; 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::Save(G4MolecularConfiguration* molconf)
{
  G4int moleculeID = molconf->GetMoleculeID();
  const G4String& moleculeName = molconf->GetFormatedName();
  G4TrackList* trackList = G4ITTrackHolder::Instance()->GetMainList(moleculeID);

  if(trackList == 0) return;

  G4TrackList::iterator it  = trackList->begin();
  G4TrackList::iterator end = trackList->end();

  for (; it != end; ++it)
  {
    G4Track* track = *it;
    SaveMoleculeInfo(track, moleculeID, moleculeName);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::SaveMoleculeInfo(
                  G4Track* track, G4int molID, const G4String& /*moleculeName*/)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if(!analysisManager->IsActive()) {return; }

  const G4ThreeVector& position = track->GetPosition();
    //H_{3}O^{1} ID = 0
    //OH^{-1}    ID = 1
    //OH^{0}     ID = 2
    //e_{aq}^{1} ID = 3  
    //H0         ID = 4
    //H_{2}^{0}  ID = 5
    //H2O2       ID = 6
    //H_{2}O^{0} or H_{2}O^{1} ID=7-17
         
    G4double xp    = position.x();
    G4double yp    = position.y();
    G4double zp    = position.z();
    G4double R     = std::sqrt(xp*xp+yp*yp+zp*zp)  /CLHEP::nm;
    //G4double RNP   = fpDetector->GetNPRadius()/CLHEP::nm;
    
    G4int offset = 10;

    if(molID<7){
      analysisManager->FillH1(molID+offset, R);
    }else{
      analysisManager->FillH1(17, R);
    }

    G4Scheduler::Instance()->Stop();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserPreTimeStepAction()
{
  // Loop over defined molecules
  G4ConfigurationIterator it 
                = G4MoleculeTable::Instance()->GetConfigurationIterator();

  G4double time = G4Scheduler::Instance()->GetGlobalTime();
  
  if(time == 1.*picosecond){
    while(it())
    {
      G4MolecularConfiguration* molconf = it.value();
      Save(molconf);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserReactionAction(
                const G4Track&,
                const G4Track&,
                const std::vector<G4Track*>* /*products*/)
{
}
