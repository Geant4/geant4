#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include"G24DetectorConstruction.hpp"
#include"G24ActionInitialization.hpp"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include"G24PhysicsList.hpp"
#include "VG01SteppingVerboseWithDir.hh"
#include "G4SystemOfUnits.hh"
#include <stdlib.h>
#include "G4PhysicalConstants.hh"
#include "G4MuonMinus.hh"
#include<iostream>

int main(int argc, char *argv[])
{
    G4double minEpsilon = 1.0e-9;
    G4double beta = 0.1; //default
    G4double step_fraction = 1.0;
    G4double mass_muon_c2 = 105.6583755*CLHEP::MeV;
    G4double mass_muon_kg =((mass_muon_c2)/c_squared)/kg;//G4MuonMinus::MuonMinusDefinition()->GetPDGMass();
    
    G4double bval_si= 1.0;//CLHEP::electron_mass_c2/mass_muon_c2; //tesla
     
    G4VSteppingVerbose::SetInstance( new VG01SteppingVerboseWithDir); 
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    char *ptr;
   
    if(argc>1) {
       minEpsilon= strtod(argv[1], nullptr);
       beta = strtod(argv[2], nullptr);  
       step_fraction = strtod(argv[3], nullptr);
   }


    G24DetectorConstruction* detectorConstructor = new G24DetectorConstruction;
    detectorConstructor-> minEpsilon = minEpsilon;
    detectorConstructor-> step_fraction = step_fraction;
    constexpr G4double mass_e_si = (electron_mass_c2/c_squared)/CLHEP::kg;;  //9.108e-31;
    G4cout<<"Mass of Muon "<<mass_muon_kg<<" kg "<<G4endl; 
    
     G4double velocity_si = beta*CLHEP::c_light/CLHEP::m*CLHEP::second;
   constexpr  G4double charge_si = CLHEP::e_SI;
   G4double gamma = 1/std::sqrt(1 - sqr(beta));
   
   G4double radius_si = (gamma*mass_muon_kg*velocity_si)/(charge_si*bval_si);
   G4cout<<"Beta: "<<beta<<G4endl;
   G4cout<<"Radius: "<<radius_si<<" m "<<G4endl;
   G4cout<<"step_fraction "<<step_fraction<<G4endl;
   G4cout<<"Speed "<<velocity_si<<G4endl;
   const G4double circum_elec =2*CLHEP::pi*radius_si*m;

     runManager->SetUserInitialization(detectorConstructor);
     detectorConstructor->circumference = circum_elec;
     detectorConstructor->bz_si = bval_si;
     //detectorConstructor->energy = 0.5*mass_e*velocity*velocity;
     
   

     //runManager->SetUserInitialization(new G24PhysicsList);
     G4VModularPhysicsList* physicsList = new FTFP_BERT;
     physicsList->RegisterPhysics(new G4StepLimiterPhysics());
     runManager->SetUserInitialization(physicsList);

    auto actionini = new G24ActionInitialization;
    actionini->fbeta = beta;
    actionini->mass_c2 = mass_muon_c2;
     runManager->SetUserInitialization(actionini);
      
    // auto visManager = new G4VisExecutive;
     //visManager->Initialize();

     runManager->Initialize();
     G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/run/verbose 1");
    UI->ApplyCommand("/event/verbose 1");
    UI->ApplyCommand("/tracking/verbose 1");

    // start a run
    int numberOfEvent = 1;
    runManager->BeamOn(numberOfEvent);

    // job termination
   // delete visManager;
    delete runManager;
    
    return 0;
    
}