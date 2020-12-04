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

#include "G4Timer.hh"
#include "GammaKnifeRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "GammaKnifeDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4ScoringManager.hh"
#include <ctime>
#include "G4SystemOfUnits.hh"


#include "GammaKnifeRunAction.hh"

////////////////////////////////////////////////////////////////////////

GammaKnifeRunAction::GammaKnifeRunAction(G4bool /*isOnMaster*/)
{
  sum_real = 0;
  sum_user = 0;
  timer = new G4Timer;
}

////////////////////////////////////////////////////////////////////////

GammaKnifeRunAction::~GammaKnifeRunAction()
{ 
  delete timer;
}

////////////////////////////////////////////////////////////////////////

void GammaKnifeRunAction::BeginOfRunAction(const G4Run* aRun)
{ 	
   G4RunManager::GetRunManager()->SetRandomNumberStore(true);
   G4cout << "Run " << aRun -> GetRunID() << " starts ..." << G4endl;

   timer->Start();  
}

////////////////////////////////////////////////////////////////////////

void GammaKnifeRunAction::EndOfRunAction(const G4Run* aRun)
{

  timer->Stop();
 
  G4double time_real = timer -> GetRealElapsed();
  G4double time_user = timer -> GetUserElapsed();
  
  sum_real = sum_real + time_real;
  sum_user = sum_user + time_user;
  G4cout << "\t  User TOT = " << sum_user << 
  "\t  Real TOT = "  << sum_real <<  G4endl;
  G4cout << " Summary of Run " << aRun -> GetRunID() <<" :"<< G4endl;

}



