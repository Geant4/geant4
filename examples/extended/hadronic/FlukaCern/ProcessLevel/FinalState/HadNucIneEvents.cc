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
///  \file HadNucIneEvents.cc
///  \brief Main program,
///         hadronic/FlukaCern/ProcessLevel/FinalState example.
//
//  Author: A. Ribbon, 8 November 2020
//  Modified: G. Hugo, 8 December 2022
//
//------------------------------------------------------------------------
//
//     HadNucIneEvents
//
/// This program is an adaptation of Hadr09 example.
/// It offers all Hadr09 features, and adds the possibility of 
/// accessing hadron-nucleus inelastic interactions final states from FLUKA.
///
/// With respect to the Hadr09 example, 
/// the program also adds the possibility of plotting the final state:
/// all encountered secondaries spectra are automatically plotted,
/// as well as the residual nuclei distributions.
/// All plots (created via the G4 analysis manager) can be dumped 
/// to any of the usually supported formats (e.g. ROOT format), 
/// but also in a Flair-compatible format.
///
/// The final states (i.e. secondary particles) produced by
/// hadron-nuclear inelastic collisions are handled by HadronicGenerator.
///
/// The use of the class Hadronic Generator is very simple:
/// the constructor needs to be invoked only once - specifying the name
/// of the "physics case" to consider ("CFLUKAHI" will be
/// considered as default if the name is not specified) - and then one
/// method needs to be called at each collision, specifying the type of
/// collision (hadron, energy, direction, material) to be simulated.
/// The class HadronicGenerator is expected to work also in a
/// multi-threaded environment with "external" threads (i.e. threads
/// that are not necessarily managed by Geant4 run-manager):
/// each thread should have its own instance of the class.
///
/// See the string "***LOOKHERE***" below for the setting of parameters
/// of this example: the "physics case", the set of possibilities from
/// which to sample the projectile
/// a list of hadrons is possible from which to sample at each collision),
/// the kinetic energy of the projectile (which can be sampled within
/// an interval), whether the direction of the projectile is fixed or
/// sampled at each collision, the target material (a list of materials
/// is possible, from which the target material can be sampled at each
/// collision, and then from this target material, the target nucleus
/// will be chosen randomly by Geant4 itself), and whether to print out
/// some information or not and how frequently.
/// Once a well-defined type of hadron-nucleus
/// inelastic collision has been chosen, the method
///   HadronicGenerator::GenerateInteraction
/// returns the secondaries produced by that interaction (in the form
/// of a G4VParticleChange object).
///
/// Here by default, an already well-defined type of hadron-nucleus
/// inelastic collision is selected 
/// (specific hadron, at a given kinetic energy and direction, 
/// on a specific material). 
/// The initial random seed is not set randomly, 
/// so that results are reproducible from one simulation to the next.
///
/// Use:  build/HadNucIneEvents
//
//------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <chrono>
#include <iomanip>

#include "globals.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VParticleChange.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "HadronicGenerator.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "FinalStateHistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int main(G4int argc, char** argv) {
  
  G4cout << "=== Test of the HadronicGenerator ===" << G4endl;

  // See the HadronicGenerator class for the possibilities and meaning of the "physics cases".
  // ( In short, it is the name of the Geant4 hadronic model used for the simulation of
  //   the collision, with the possibility of having a transition between two models in
  //   a given energy interval, as in physics lists. )

  //***LOOKHERE***  PHYSICS CASE
  G4String namePhysics = "CFLUKAHI";
  //const G4String namePhysics = "FTFP_BERT";       
  //const G4String namePhysics = "FTFP_BERT_ATL";
  //const G4String namePhysics = "QGSP_BERT";
  //const G4String namePhysics = "QGSP_BIC";
  //const G4String namePhysics = "FTFP_INCLXX";
  //const G4String namePhysics = "FTFP";
  //const G4String namePhysics = "QGSP";
  //const G4String namePhysics = "BERT";
  //const G4String namePhysics = "BIC";
  //const G4String namePhysics = "IonBIC";
  //const G4String namePhysics = "INCL";
  
  // The kinetic energy of the projectile will be sampled randomly, with flat probability
  // in the interval [minEnergy, maxEnergy].
  G4double minEnergy = 7.*CLHEP::TeV;  //***LOOKHERE***  HADRON PROJECTILE MIN Ekin
  G4double maxEnergy = 7.*CLHEP::TeV;  //***LOOKHERE***  HADRON PROJECTILE MAX Ekin
  
  G4int numCollisions = 100000;  //***LOOKHERE***  NUMBER OF COLLISIONS
  //const G4int numCollisions = 100;    // DEBUG

  // IMPORTANT - TESTING ONLY: 
  // OVERWRITES DEFAULT PHYSICS CASE AND NUMBER OF EVENTS
  std::vector<G4String> args(argv, argv + argc);
  if (args.size() == 2 && args[1] == "--test") {
    namePhysics = G4String("FTFP_BERT");
    numCollisions = 10;
  }

  // Enable or disable the print out of this program: if enabled, the number of secondaries
  // produced in each collisions is printed out; moreover, once every "printingGap"
  // collisions, the list of secondaries is printed out.
  const G4bool isPrintingEnabled = true;           //***LOOKHERE***  PRINT OUT ON/OFF
  const G4int  printingGap = 100;                  //***LOOKHERE***  GAP IN PRINTING
  
  // Vector of Geant4 names of hadron projectiles: one of this will be sampled randomly
  // (with uniform probability) for each collision, when the projectile is not an ion.
  // Note: comment out the corresponding line in order to exclude a particle.
  std::vector< G4String > vecProjectiles;  //***LOOKHERE***  POSSIBLE HADRON PROJECTILES
  //vecProjectiles.push_back( "pi-" );
  //Note: vecProjectiles.push_back( "pi0" );               // Excluded because too short-lived
  //vecProjectiles.push_back( "pi+" );
  //vecProjectiles.push_back( "kaon-" );
  //vecProjectiles.push_back( "kaon+" );
  //vecProjectiles.push_back( "kaon0L" );
  //vecProjectiles.push_back( "kaon0S" );
  //Note: vecProjectiles.push_back( "eta" );               // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "eta_prime" );         // Excluded because too short-lived
  vecProjectiles.push_back( "proton" );
  //vecProjectiles.push_back( "neutron" );
  //vecProjectiles.push_back( "deuteron" );
  //vecProjectiles.push_back( "triton" );
  //vecProjectiles.push_back( "He3" );
  //vecProjectiles.push_back( "alpha" );
  //vecProjectiles.push_back( "lambda" );
  //vecProjectiles.push_back( "sigma-" );
  //Note: vecProjectiles.push_back( "sigma0" );            // Excluded because too short-lived
  //vecProjectiles.push_back( "sigma+" );
  //vecProjectiles.push_back( "xi-" );
  //vecProjectiles.push_back( "xi0" );
  //vecProjectiles.push_back( "omega-" );
  //vecProjectiles.push_back( "anti_proton" );
  //vecProjectiles.push_back( "anti_neutron" );
  //vecProjectiles.push_back( "anti_lambda" );
  //vecProjectiles.push_back( "anti_sigma-" );
  //Note: vecProjectiles.push_back( "anti_sigma0" );       // Excluded because too short-lived
  //vecProjectiles.push_back( "anti_sigma+" );
  //vecProjectiles.push_back( "anti_xi-" );
  //vecProjectiles.push_back( "anti_xi0" );
  //vecProjectiles.push_back( "anti_omega-" );
  //vecProjectiles.push_back( "anti_deuteron" );
  //vecProjectiles.push_back( "anti_triton" );
  //vecProjectiles.push_back( "anti_He3" );
  //vecProjectiles.push_back( "anti_alpha" );
  // Charm and bottom hadrons
  //vecProjectiles.push_back( "D+" );
  //vecProjectiles.push_back( "D-" );
  //vecProjectiles.push_back( "D0" );
  //vecProjectiles.push_back( "anti_D0" );
  //vecProjectiles.push_back( "Ds+" );
  //vecProjectiles.push_back( "Ds-" );
  //Note: vecProjectiles.push_back( "etac" );              // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "J/psi" );             // Excluded because too short-lived
  //vecProjectiles.push_back( "B+" );
  //vecProjectiles.push_back( "B-" );
  //vecProjectiles.push_back( "B0" );
  //vecProjectiles.push_back( "anti_B0" );
  //vecProjectiles.push_back( "Bs0" );
  //vecProjectiles.push_back( "anti_Bs0" );
  //vecProjectiles.push_back( "Bc+" );
  //vecProjectiles.push_back( "Bc-" );
  //Note: vecProjectiles.push_back( "Upsilon" );           // Excluded because too short-lived
  //vecProjectiles.push_back( "lambda_c+" );
  //vecProjectiles.push_back( "anti_lambda_c+" );
  //Note: vecProjectiles.push_back( "sigma_c+" );          // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "anti_sigma_c+" );     // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "sigma_c0" );          // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "anti_sigma_c0" );     // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "sigma_c++" );         // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "anti_sigma_c++" );    // Excluded because too short-lived
  //vecProjectiles.push_back( "xi_c+" );
  //vecProjectiles.push_back( "anti_xi_c+" );
  //vecProjectiles.push_back( "xi_c0" );
  //vecProjectiles.push_back( "anti_xi_c0" );
  //vecProjectiles.push_back( "omega_c0" );
  //vecProjectiles.push_back( "anti_omega_c0" );
  //vecProjectiles.push_back( "lambda_b" );
  //vecProjectiles.push_back( "anti_lambda_b" );
  //Note: vecProjectiles.push_back( "sigma_b+" );          // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "anti_sigma_b+" );     // Excluded because too short-lived  
  //Note: vecProjectiles.push_back( "sigma_b0" );          // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "sigma_b0" );          // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "sigma_b-" );          // Excluded because too short-lived
  //Note: vecProjectiles.push_back( "anti_sigma_b-" );     // Excluded because too short-lived  
  //vecProjectiles.push_back( "xi_b0" );
  //vecProjectiles.push_back( "anti_xi_b0" );
  //vecProjectiles.push_back( "xi_b-" );
  //vecProjectiles.push_back( "anti_xi_b-" );
  //vecProjectiles.push_back( "omega_b-" );
  //vecProjectiles.push_back( "anti_omega_b-" );

  G4ParticleDefinition* projectileNucleus = nullptr;
  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager( new G4ProcessManager( gion ) );
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ions = partTable->GetIonTable();
  partTable->SetReadiness();
  ions->CreateAllIon();
  ions->CreateAllIsomer();
  
  //***LOOKHERE***  HADRON (false) OR ION (true) PROJECTILE ?
  const G4bool isProjectileIon = false;
  if ( isProjectileIon ) {
    minEnergy = 40.0*13.0*CLHEP::GeV;    //***LOOKHERE***  ION PROJECTILE MIN Ekin
    maxEnergy = 40.0*13.0*CLHEP::GeV;    //***LOOKHERE***  ION PROJECTILE MAX Ekin
    G4int ionZ = 18, ionA = 40;          //***LOOKHERE***  ION PROJECTILE (Z, A)
    projectileNucleus = partTable->GetIonTable()->GetIon( ionZ, ionA, 0.0 );
  }

  // Vector of Geant4 NIST names of materials: one of this will be sampled randomly
  // (with uniform probability) for each collision and used as target material.
  // Note: comment out the corresponding line in order to exclude a material;
  //       or, vice versa, add a new line to extend the list with another material.
  std::vector< G4String > vecMaterials;  //***LOOKHERE*** : NIST TARGET MATERIALS
  //vecMaterials.push_back( "G4_H" );
  //vecMaterials.push_back( "G4_He" );
  //vecMaterials.push_back( "G4_Be" );
  vecMaterials.push_back( "G4_C" );
  //vecMaterials.push_back( "G4_Al" );
  //vecMaterials.push_back( "G4_Si" );
  //vecMaterials.push_back( "G4_Sc" );
  //vecMaterials.push_back( "G4_Ar" );
  //vecMaterials.push_back( "G4_Fe" );
  //vecMaterials.push_back( "G4_Cu" );
  //vecMaterials.push_back( "G4_W" );
  //vecMaterials.push_back( "G4_Pb" );
  
  const G4int numProjectiles = vecProjectiles.size();
  const G4int numMaterials = vecMaterials.size();

  G4cout << G4endl
         << "=================  Configuration ==================" << G4endl
         << "Model: " << namePhysics << G4endl
         << "Ekin: [ " << minEnergy/CLHEP::GeV << " , " << maxEnergy/CLHEP::GeV
         << " ] GeV" << G4endl
         << "Number of collisions:  " << numCollisions << G4endl
         << "Number of hadron projectiles: " << numProjectiles << G4endl
         << "Number of materials:   " << numMaterials   << G4endl
         << "IsIonProjectile: " << ( projectileNucleus != nullptr ? "true \t" : "false" )
         << ( projectileNucleus != nullptr ? projectileNucleus->GetParticleName() : "") << G4endl
         << "===================================================" << G4endl
         << G4endl;
  
  CLHEP::Ranlux64Engine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  //***LOOKHERE***  RANDOM ENGINE START SEED
  //G4int seed = time( NULL ); 
  //CLHEP::HepRandom::setTheSeed( seed ); 
  //G4cout << G4endl << " Initial seed = " << seed << G4endl << G4endl; 

  // Set up histo manager.
  auto histoManager = FinalStateHistoManager();
  histoManager.Book();
  
  // Instanciate the HadronicGenerator providing the name of the "physics case"
  HadronicGenerator* theHadronicGenerator = new HadronicGenerator( namePhysics );
  //****************************************************************************
  
  if ( theHadronicGenerator == nullptr ) {
    G4cerr << "ERROR: theHadronicGenerator is NULL !" << G4endl;
    return 1;
  } else if ( ! theHadronicGenerator->IsPhysicsCaseSupported() ) {
    G4cerr << "ERROR: this physics case is NOT supported !" << G4endl;
    return 2;
  }
  

  // Start timing
  auto start = std::chrono::high_resolution_clock::now();

  // Loop over the collisions
  G4double rnd1, rnd2, rnd3, rnd4, rnd5, rnd6, normalization, projectileEnergy;
  G4VParticleChange* aChange = nullptr;
  for ( G4int i = 0; i < numCollisions; ++i ) {
    histoManager.BeginOfEvent();

    // Draw some random numbers to select the hadron-nucleus interaction:
    // projectile hadron, projectile kinetic energy, projectile direction, and target material.
    rnd1 = CLHEP::HepRandom::getTheEngine()->flat(); 
    rnd2 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd3 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd4 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd5 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd6 = CLHEP::HepRandom::getTheEngine()->flat();
    // Sample the projectile kinetic energy
    projectileEnergy = minEnergy + rnd1*( maxEnergy - minEnergy );
    if ( projectileEnergy <= 0.0 ) projectileEnergy = minEnergy; 
    // Sample the projectile direction
    normalization = 1.0 / std::sqrt( rnd2*rnd2 + rnd3*rnd3 + rnd4*rnd4 );
    //***LOOKHERE***  IF true THEN SMEAR DIRECTION 
    const G4bool isOnSmearingDirection = false ;
    //***LOOKHERE***  ELSE USE THIS FIXED DIRECTION
    G4ThreeVector aDirection = G4ThreeVector( 0.0, 0.0, 1.0 );
    if ( isOnSmearingDirection ) {
      aDirection = G4ThreeVector( normalization*rnd2, normalization*rnd3, normalization*rnd4 );
    } 
    // Sample the projectile hadron from the vector vecProjectiles
    G4int index_projectile = std::trunc( rnd5*numProjectiles );
    G4String nameProjectile = vecProjectiles[ index_projectile ];
    G4ParticleDefinition* projectile = partTable->FindParticle( nameProjectile );
    if ( projectileNucleus ) {
      nameProjectile = projectileNucleus->GetParticleName();
      projectile = projectileNucleus;
    }
    // Sample the target material from the vector vecMaterials
    // (Note: the target nucleus will be sampled by Geant4)
    G4int index_material = std::trunc( rnd6*numMaterials );
    G4String nameMaterial = vecMaterials[ index_material ];
    G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial( nameMaterial );
    if ( material == nullptr ) {
      G4cerr << "ERROR: Material " << nameMaterial << " is not found !" << G4endl;
      return 3;
    }
    if ( isPrintingEnabled ) {
      G4cout << "\t Collision " << i << " ; projectile=" << nameProjectile;
      if ( projectileNucleus ) {
        G4cout << " ; Ekin[MeV]/nucleon=" << projectileEnergy / 
          static_cast< G4double >( std::abs( projectileNucleus->GetBaryonNumber() ) );
      } else {
        G4cout << " ; Ekin[MeV]=" << projectileEnergy; 
      }
      G4cout << " ; direction=" << aDirection << " ; material=" << nameMaterial;
    }
    
    // Call here the "hadronic generator" to get the secondaries produced by the hadronic collision
    aChange = theHadronicGenerator->GenerateInteraction( projectile, projectileEnergy,
    /* ********************************************** */ aDirection, material );
    
    G4int nsec = aChange ? aChange->GetNumberOfSecondaries() : 0;
    G4bool isPrintingOfSecondariesEnabled = false;
    if ( isPrintingEnabled ) {
      G4cout << G4endl << "\t --> #secondaries=" << nsec 
             << " ; impactParameter[fm]=" 
             << theHadronicGenerator->GetImpactParameter() / fermi
             << " ; #projectileSpectatorNucleons=" 
             << theHadronicGenerator->GetNumberOfProjectileSpectatorNucleons()
             << " ; #targetSpectatorNucleons=" 
             << theHadronicGenerator->GetNumberOfTargetSpectatorNucleons()
             << " ; #NNcollisions=" 
             << theHadronicGenerator->GetNumberOfNNcollisions() 
             << G4endl;
      if ( i % printingGap == 0 ) {
        isPrintingOfSecondariesEnabled = true;
        G4cout << "\t \t List of produced secondaries: " << G4endl;
      }
    }
    // Loop over produced secondaries and eventually print out some information.
    for ( G4int j = 0; j < nsec; ++j ) {
      const G4DynamicParticle* sec = aChange->GetSecondary(j)->GetDynamicParticle();
      if ( isPrintingOfSecondariesEnabled ) {
        G4cout << "\t \t \t j=" << j << "\t" << sec->GetDefinition()->GetParticleName()
               << "\t p=" << sec->Get4Momentum() << " MeV" << G4endl;
      }

      // Store each secondary.
      histoManager.ScoreSecondary(sec);

      delete aChange->GetSecondary(j);
    }
    if ( aChange ) aChange->Clear();
    histoManager.EndOfEvent();
  }

  histoManager.EndOfRun();


  G4cout << G4endl << " Final random number = " << CLHEP::HepRandom::getTheEngine()->flat() 
         << G4endl;

  const auto stop = std::chrono::high_resolution_clock::now();
  const auto diff = stop - start;
  const auto time = static_cast<G4double>(
      std::chrono::duration_cast<std::chrono::microseconds>(diff).count()) / 1e6;
  G4cout << G4endl;
  G4cout << "Processed " << numCollisions << " events (collisions) in " 
         << std::scientific << time << " seconds."
         << " Average: " << std::defaultfloat << (time * 1E3 / numCollisions) << " ms / event." 
         << G4endl;
  G4cout << G4endl;

  G4cout << "=== End of test ===" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
