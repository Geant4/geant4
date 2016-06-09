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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISAEventAction class                                            *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************


#include "LISAEventAction.hh"

#include "LISAPrimaryGeneratorAction.hh"
#include "LISASteppingAction.hh"

#ifdef G4ANALYSIS_USE
#include "LISAAnalysisManager.hh"
#endif

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ParticleDefinition.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include <fstream>
#include <iomanip>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISAEventAction::LISAEventAction(LISAPrimaryGeneratorAction* generatorAction,
  LISASteppingAction* steppingAction)
  : genAction(generatorAction),stepAction(steppingAction) {

  energy_pri    = 0;
  charge_in[0]  = 0;
  charge_in[1]  = 0;
  charge_out[0] = 0;
  charge_out[1] = 0;
  charge_tot[0] = 0;
  charge_tot[1] = 0;
  seeds         = NULL;
  filename      = G4String();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISAEventAction::~LISAEventAction() {;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISAEventAction::BeginOfEventAction(const G4Event* evt) {

  // reset event charges
  stepAction->Discharge();

  // grab event seeds
  seeds = genAction->GetEventSeeds();
  //    G4cout << " 1st seed: " << *seeds << G4endl;;
  //    G4cout << " 2nd seed: " << *(seeds+1) << G4endl;
  //    HepRandom::showEngineStatus();

  // grab energy of primary
  energy_pri = genAction->GetEnergyPrimary();

  event_id = evt->GetEventID();
  if(event_id%1000==0 || (event_id%100==0 && event_id<1000))
    G4cout << "\n---> Begin of event: " << event_id << G4endl;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISAEventAction::EndOfEventAction(const G4Event* evt) {


  // Look at TM0 and TM1
  for(int tm=0; tm<=1; ++tm) {

    charge_in[tm]  = stepAction->GetChargeIn(tm);
    charge_out[tm] = stepAction->GetChargeOut(tm);
    charge_tot[tm] = charge_in[tm] - charge_out[tm];

    // save if changed
    if(charge_tot[tm]) {

      // printout
      G4cout << " *** Evt: " << event_id
	     << "\tTM: "     << tm
	     << "\tE: "      << energy_pri
	     << "\tQnet: "   << charge_tot[tm]
	     << "\tQin: "    << charge_in[tm] 
	     << "\tQout: "   << charge_out[tm] 
	     << G4endl;
    
      // append to file
      std::ofstream chargefile(filename,std::ios::app);
      chargefile << event_id       << "\t" 
		 << tm             << "\t"
		 << energy_pri     << "\t" 
		 << charge_tot[tm] << "\t"
		 << charge_in[tm]  << "\t"
		 << charge_out[tm] << "\t"
		 << *seeds         << "\t"
		 << *(seeds+1)     << "\t"
		 << G4endl;
      
      // Write event information to HBOOK NTuple if using analysis
#ifdef G4ANALYSIS_USE
      LISAAnalysisManager* analysis = LISAAnalysisManager::getInstance();
      analysis->analyseRun(event_id, 
			   tm,
			   energy_pri,
			   charge_tot[tm],
			   charge_in[tm],
			   charge_out[tm],
			   *seeds, 
			   *(seeds+1));
#endif

    } // if(charge_tot[tm])

  } // for(int tm=0; tm<=1; ++tm)



  // draw trajectories
  // colours by particle type
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {

    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) 
      n_trajectories = trajectoryContainer->entries();

    for (G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* trj = 0;
      trj = dynamic_cast<G4Trajectory*>((*trajectoryContainer)[i]);

      G4String name = trj->GetParticleDefinition()->GetParticleName();
      G4Colour colour = G4Colour(1.0, 1.0, 1.0);
      G4Polyline pPolyline;

      if (trj) {
	
	for (int i=0; i<trj->GetPointEntries(); i++)
	  pPolyline.push_back( trj->GetPoint(i)->GetPosition() );
	
	     if(name=="proton")   colour=G4Colour(0.0, 0.0, 1.0);
	else if(name=="alpha")    colour=G4Colour(0.0, 0.0, 0.5);
	else if(name=="He3")      colour=G4Colour(1.0, 0.0, 1.0);
	else if(name=="pi+")      colour=G4Colour(0.5, 0.5, 0.5);
	else if(name=="pi-")      colour=G4Colour(0.5, 0.5, 0.5);
	else if(name=="e+")       colour=G4Colour(0.0, 1.0, 1.0);
	else if(name=="e-")       colour=G4Colour(1.0, 0.0, 0.0);
	else if(name=="gamma")    colour=G4Colour(0.0, 1.0, 0.0);
	else if(name=="neutron")  colour=G4Colour(1.0, 1.0, 0.0);
	else                      colour=G4Colour(1.0, 1.0, 1.0);

	G4VisAttributes attribs(colour);
	pPolyline.SetVisAttributes(attribs);
	if(pVVisManager) pVVisManager->Draw(pPolyline);
	
      } else
	G4cerr << "EndOfEventAction: Failed to dynamic cast to G4Trajectory!" 
	       << G4endl;
    }
  }


}

