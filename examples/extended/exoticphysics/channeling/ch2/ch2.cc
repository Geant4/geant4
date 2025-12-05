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
/// \file ch2.cc
/// \brief Main program of the channeling/ch2 example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4FastSimulationPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "G4Timer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Get current time
  G4Timer* theTimer = new G4Timer();
  theTimer->Start();
  
  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(0.);
  
  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }
  
  // Construct the default run manager
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new DetectorConstruction());

  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;

  // -- Create helper tool, used to activate the fast simulation:
  G4FastSimulationPhysics* fastSimulationPhysics = new G4FastSimulationPhysics();
  fastSimulationPhysics->BeVerbose();
  // -- activation of fast simulation for particles having fast simulation models
  // -- attached in the mass geometry:
  // you may add any charged particles here
  // CAUTION: for the particles other then e+- you would likely want
  // to switch off the radiation
  fastSimulationPhysics->ActivateFastSimulation("e-");
  fastSimulationPhysics->ActivateFastSimulation("e+");
  fastSimulationPhysics->ActivateFastSimulation("proton");
  fastSimulationPhysics->ActivateFastSimulation("anti_proton");
  fastSimulationPhysics->ActivateFastSimulation("mu+");
  fastSimulationPhysics->ActivateFastSimulation("mu-");
  fastSimulationPhysics->ActivateFastSimulation("pi+");
  fastSimulationPhysics->ActivateFastSimulation("pi-");
  fastSimulationPhysics->ActivateFastSimulation("GenericIon");
  //fastSimulationPhysics->ActivateFastSimulation("your_particle");
  // you may activate this model for any charged particle
  // a neutral particle will not enter the model

  // -- Attach the fast simulation physics constructor to the physics list:
  physicsList->RegisterPhysics( fastSimulationPhysics );
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization());

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  G4int nofEventsTot = runManager->GetNumberOfEventsToBeProcessed();

  delete visManager;
  delete runManager;

  //final output for the spectrum
  std::vector<G4double> photonEnergyInSpectrum;
  std::vector<G4double> spectrum;
  G4bool spectrumWrite = false;

  for (G4int ii = 0; ; ++ii)
  {
      std::string filename = "Spectrum_"+std::to_string(ii)+".dat";
      std::ifstream fileN1(filename);

      //if no file
      if (!fileN1) break;
      spectrumWrite = true; //if a temporary file exist => final output will be written

      G4int jj = 0;
      G4double eph, spec=0.;
      while (fileN1 >> eph >> spec)
      {
          if (ii==0)
          {
              photonEnergyInSpectrum.push_back(eph);//the same for all the files
              spectrum.push_back(spec); //the same data size for all the files
          }
          else{spectrum[jj++] += spec;} // spectrum accumulation
      }
      fileN1.close();

      std::remove(filename.c_str());
  }

  std::remove("Spectrum.dat");//delete a previous file if existed
  if(spectrumWrite)// if radiation model = false => no spectrum
  {
      std::ofstream file1;
      file1.open("Spectrum.dat");
      //CAUTION: spectrum is normalized onto a probability of radiation W, IT IS NOT A PDF
      file1 << "# E_photon [MeV]" <<
          "  dW_rad/dE_photon [MeV^-1]; W_rad - radiation probability" << G4endl;
      for(std::size_t i = 0; i<spectrum.size(); i++)
      {file1 << photonEnergyInSpectrum[i] << " " << spectrum[i]/nofEventsTot << G4endl;}
      file1.close();
  }

  theTimer->Stop();
  G4cout << "Execution terminated" << G4endl;
  G4cout << (*theTimer) << G4endl;
  delete theTimer;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
