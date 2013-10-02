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
// $Id: B1ConRun.cc 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file B1ConRun.cc
/// \brief Implementation of the B1ConRun class

#include "B1ConRun.hh"
#include "B1EventInformation.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ConRun::B1ConRun()
: G4Run(),
  fEdepRun(0.), fEdep2Run(0.)
{ fEdepEventVector.clear(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ConRun::~B1ConRun()
{ ; } 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ConRun::RecordEvent(const G4Event* event)
{  
  B1EventInformation* evInfo
    = static_cast<B1EventInformation*>(event->GetUserInformation());
  G4double EdepEvent = evInfo->GetEdepEvent();
  fEdepRun  += EdepEvent;
  fEdep2Run += EdepEvent*EdepEvent;
  fEdepEventVector.push_back( evInfo->GetEdepEvent() );

  G4Run::RecordEvent(event);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ConRun::Merge(const G4Run* aRun)
{
  const B1ConRun* localRun = static_cast<const B1ConRun*>(aRun);
  //const B1ConRun* localRunn = static_cast<const B1ConRun*>(aRun);
  //B1ConRun* localRun = const_cast<B1ConRun*>(localRunn);

  fEdepRun  += localRun->fEdepRun;
  fEdep2Run += localRun->fEdep2Run;

  for ( size_t i = 0 ; i != localRun->fEdepEventVector.size() ; i++ ) {
     fEdepEventVector.push_back( localRun->fEdepEventVector[i] );
  }
  //for ( std::vector<G4double>::iterator 
  //      it = localRun->fEdepEventVector.begin(); it != localRun->fEdepEventVector.end(); it++ ) { 
  //   fvEdepEventVector.push_back( *it );
  //}
  G4Run::Merge(aRun); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
