//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: HadrontherapyRunAction.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include <iomanip.h>
#include "Randomize.hh"

#ifdef G4ANALYSIS_USE
#include "HadrontherapyAnalysisManager.hh"
#endif

// ---------------------------------------------------------------
HadrontherapyRunAction::HadrontherapyRunAction()
:NbOfLayer(20000)  
{
}

// ---------------------------------------------------------------
HadrontherapyRunAction::~HadrontherapyRunAction()
{
  //  delete runMessenger;
}

// ---------------------------------------------------------------
void HadrontherapyRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  for (G4int i=0; i < 50000; i++) 
    {
      energy[i] = 0.0; 
    };
  
//   ---------- Analysis preliminary implentation ------------------
// for the registration of the Bragg curve file in a .xml file

#ifdef G4ANALYSIS_USE
HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::getInstance();
analysis->book();
#endif   

G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
// save Rndm status
G4RunManager::GetRunManager()->SetRandomNumberStore(true);
HepRandom::showEngineStatus();
}

// ---------------------------------------------------------------------
void HadrontherapyRunAction::EndOfRunAction(const G4Run*) 
{  
  // WRITE ASCII FILE FOR THE REGISTRATION OF THE BRAGG PEAK
  // The a two column file (BraggPeak.out) is registered. 
  //The first column represents
  // the ionization chamber position (in mm of water); the second
  // is the energy deposited for each position of the chamber

  for (G4int i = 0; i < 20000; i++) {
    
G4double depth = 0;
depth = i*0.002; //the number represents the thickness of the ionization chamber (in mm ofwater)
    
 std::ofstream pmtfile("BraggPeak.out", std::ios::app);
 if(pmtfile.is_open())
   
   {
     pmtfile << depth  << '\t' << energy[i]  << '\t'  << G4endl; 
   } 
  }
  
  
#ifdef G4ANALYSIS_USE  
  HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::getInstance();

  //histogram fill

  analysis -> Energy_Event( slice, energy[ slice ]);

  for (slice = 0; slice < 200; slice++)  
  // It is necessary to set here the Number of slice
    {
      //G4cout << "%%%%%%%" << slice << "%%%%%" << 
      //energy[ slice ] << "%%%%%" << G4endl;
      
      //n-tuple fill
      analysis -> Energy_Dep( slice + 0.2, energy[ slice ]);         
    }
  analysis -> finish();

#endif
  
  // show Rndm status
  HepRandom::showEngineStatus();
}

// -----------------------------------------------------------------------------
void HadrontherapyRunAction::EnergyTotSlice(G4int slice, G4double energy_dep) 
{
  energy[ slice ] = energy [ slice ] + energy_dep; 
  
  //G4cout << "energy dalla funzione EnergyTotSlice" << energy[slice]
  //<< "slice" << slice << G4endl;
 
}
/*void HadrontherapyRunAction::EnergyTotMarkus()
{
 
  for (G4int i = 0; i < 250; i++) {
      energyMarkus[i] = 0;
    }
  
  for (i = 0; i < 191; i++) {
      for (G4int j = 0; j < 10; j++) {
	  energyMarkus[ i ] = energy[ i+j ] + energyMarkus[ i ];
	}   
    }
  
}
*/












