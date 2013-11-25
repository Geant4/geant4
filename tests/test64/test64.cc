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
//------------------------------------------------------------------------
// Last update: 11-Oct-2012
// 
// Author: Alberto Ribon
//
// Description: Process-level test of Bertini and Fritiof
//              nuclear capture at rest.
//
// The program expects a single input parameter: the name of a text file
// where each not empty and not commented (#) line specifies 4 things:
//   1) model: BERT or FTFP
//   2) particle type: pi-, kaon-, sigma- (for BERT);
//                     anti_proton, anti_sigma+ (for FTFP).
//   3) target material: NIST material, e.g. G4_Cu
//   4) number of events
// If the line has fewer than 4 arguments, then the line is ignored;
// if the line has more than 4 arguments, then the extra arguments
// are ignored.
//
//------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include "Randomize.hh"
#include "Tst64DetectorConstruction.hh"
#include "Tst64PhysicsList.hh"
#include "Tst64ActionInitialization.hh"
#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"

#include "G4DecayPhysics.hh"
#include "G4ParticleTable.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChange.hh"
#include "G4TrackStatus.hh"
#include "G4DynamicParticle.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"

int main( int argc, char** argv ) {

  if ( argc != 2  ) {
    G4cout << "Error: wrong number of arguments!" << G4endl
           << "Usage: mymain filename" << G4endl
           << G4endl;
    return 1;
  }

  std::string inputFileName( argv[1] );

  std::ifstream inputFile( inputFileName.c_str() );
  if ( ! inputFile ) {
    std::cerr << " Can't open input file: " << inputFileName << std::endl;
    return 2;
  }

  // Initialization

  CLHEP::Ranlux64Engine defaultEngine( 1234567, 4 ); 
  G4Random::setTheEngine(&defaultEngine);
  G4int seed = std::time( NULL ); 
  G4Random::setTheSeed(seed);
  G4cout << " Initial seed = " << seed << G4endl; 
  //G4Random::getTheEngine()->restoreStatus( "start.rndm" );

  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization( new Tst64DetectorConstruction );
  runManager->SetUserInitialization( new Tst64PhysicsList );
  runManager->SetUserInitialization(new Tst64ActionInitialization);

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4DecayPhysics decays;
  decays.ConstructParticle();

  runManager->Initialize();

  G4NistManager* nistManager = G4NistManager::Instance();

  // Parse each line of the input file

  const std::string delims( " \t" );
  std::string inputLine = "";
  while ( getline( inputFile, inputLine ) ) { 

    G4cout << "Parsing line: " << inputLine << G4endl;
    std::string::size_type begIdx, valIdx;
    begIdx = inputLine.find_first_not_of( delims );
    valIdx = inputLine.find_first_of( "#" );
    if ( begIdx == std::string::npos || begIdx == valIdx ) continue; 
    std::string modelName, particleName, materialName, ignored;
    int numberOfEvents = 0;
    std::istringstream myStream( inputLine );
    myStream >> modelName >> particleName >> materialName >> numberOfEvents >> ignored;
    if ( particleName.size() == 0 || materialName.size() == 0 || numberOfEvents == 0 ) {
      G4cout << " Too few input arguments: line skipped!" << G4endl;
      continue;
    } 
    if ( ignored.size() > 0 ) G4cout << " Ignored extra parameters!" << G4endl;
    G4cout << " ==========  Required RUN  ============" << G4endl
           << " Model = " << modelName << G4endl
           << " Particle = " << particleName << G4endl
           << " Target = " << materialName << G4endl
           << " Nevents = " << numberOfEvents << G4endl;

    // Prepare the run

    G4Material* material = 0;
    if ( nistManager ) {
      material = nistManager->FindOrBuildMaterial( materialName );
      if ( ! material ) {
        G4cerr << "ERROR: material not recognized! " << materialName << G4endl;
        continue;
      }
    }

    G4HadronStoppingProcess* processNew = 0;
    if ( modelName == "BERT" || modelName == "bert" ) {
      processNew = new G4HadronicAbsorptionBertini();
    } else if ( modelName == "FTFP" || modelName == "ftfp" ) {
      processNew = new G4HadronicAbsorptionFritiof();
    } else {
      G4cerr << "ERROR: model not recognized! " << modelName << G4endl;
      continue;
    }

    G4ParticleDefinition* part = 0;
    if ( particleName == "mu-" ) {
      part = G4MuonMinus::MuonMinus();
    } else if ( particleName == "pi-" ) {
      part = G4PionMinus::PionMinus(); 
    } else if ( particleName == "kaon-" ) {
      part = G4KaonMinus::KaonMinus();
    } else if ( particleName == "anti_proton" ) {    
      part = G4AntiProton::AntiProton();
    } else if ( particleName == "sigma-" ) {
      part = G4SigmaMinus::SigmaMinus();
    } else if ( particleName == "anti_sigma+" ) {
      part = G4AntiSigmaPlus::AntiSigmaPlus();
    } else if ( particleName == "xi-" ) {
      part = G4XiMinus::XiMinus();
    } else if ( particleName == "omega-" ) {
      part = G4OmegaMinus::OmegaMinus();
    }
    if ( ! part ) {
      G4cerr << "ERROR: particle not recognized! " << particleName << G4endl;
      continue;
    }

    G4ThreeVector aMomentum = G4ThreeVector( 0.0, 0.0, 0.0 );
    G4ThreeVector aPosition( 0.0,  0.0, 0.0 );
    G4double aTime = 0.0;

    G4StepPoint* aPoint = new G4StepPoint();
    aPoint->SetPosition( aPosition );

    G4Navigator* pNavigator = 
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    pNavigator->LocateGlobalPointAndSetup( aPosition );
    G4TouchableHandle theTouchableHandle = pNavigator->CreateTouchableHistory();  
    aPoint->SetTouchableHandle( theTouchableHandle );

    aPoint->SetMaterial( material );

    G4cout << " \t --------- Valid requested run ---------- " << G4endl
           << " \t Name process : " << processNew->GetProcessName() << G4endl
           << " \t Particle type : " << part->GetParticleName() << G4endl
           << " \t Material : " << aPoint->GetMaterial()->GetName() << G4endl
           << " \t Number of events : " << numberOfEvents << G4endl
           << " \t ---------------------------------------- " << G4endl;

    G4DynamicParticle dParticle( part, aMomentum );
    G4Track* track = new G4Track( &dParticle, aTime, aPosition );
    G4Step step;
    step.SetTrack( track );
    step.SetPreStepPoint( aPoint ); 
    track->SetStep( &step );
    track->SetTouchableHandle( theTouchableHandle );

    G4double initialTotalEnergy = dParticle.GetKineticEnergy() / GeV;
    if ( part->GetParticleName() == "pi+"     ||
         part->GetParticleName() == "pi-"     ||
         part->GetParticleName() == "pi0"     ||
         part->GetParticleName() == "kaon+"   ||
         part->GetParticleName() == "kaon-"   ||
         part->GetParticleName() == "kaon0S"  ||
         part->GetParticleName() == "kaon0L" ) {
      initialTotalEnergy += dParticle.GetMass() / GeV;
    } else if ( part->GetParticleName() == "anti_proton" ) {
      initialTotalEnergy += 2.0 * dParticle.GetMass() / GeV;
    } else if ( part->GetParticleName() == "sigma-"  ||
                part->GetParticleName() == "xi-"     ||
                part->GetParticleName() == "omega-" ) {
      initialTotalEnergy += dParticle.GetMass() / GeV - 0.939;
    } else if ( part->GetParticleName() == "anti_sigma+" ) {
      initialTotalEnergy += dParticle.GetMass() / GeV + 0.939;
    }
    G4cout << "\t Initial total energy = " << initialTotalEnergy << " GeV" << G4endl;  

    G4VParticleChange* aChange = 0;
    G4Track* sec = 0;
    G4String secName;
    G4double secKE = 0.0;
    G4double secTE = 0.0;
    G4double secMass = 0.0;

    G4int cumulatedNumberSecondaries = 0;
    G4double cumulatedSumEnergySecondaries = 0.0;
    G4int countOverProduction = 0;
    G4int countUnderProduction = 0;
    G4double maxSumEnergySecondaries = 0.0;
    G4double minSumEnergySecondaries = 999.9; 

    G4cout << "\t Start event loop..." << G4endl;

    for ( G4int loop = 0; loop < numberOfEvents; loop++ ) {     // Event loop

      if ( loop % 1000 == 0 ) {
        G4cout << "\t \t -> Event " << loop << G4endl;
      }

      if ( loop % 100000 == 0 ) {
        std::stringstream out;
        out << "currentEvent_" << loop << ".rndm";
        std::string nameFile( out.str() );
        G4Random::getTheEngine()->saveStatus( nameFile.c_str() );
      }

      if ( processNew )  {
        aChange = processNew->AtRestDoIt( *track, step );
      }

      G4int nSec = aChange->GetNumberOfSecondaries();
      G4double sumEnergySecondaries = 0.0;

      for ( G4int j = 0; j < nSec; j++ ) {    // Produced secondary loop

        sec = aChange->GetSecondary( j );
        secName = sec->GetDefinition()->GetParticleName();
        secKE = sec->GetKineticEnergy()/GeV;
        secTE = sec->GetTotalEnergy()/GeV;
        secMass = sec->GetDefinition()->GetPDGMass()/GeV;

        if ( secName == "gamma"  || 
             secName == "pi+"    ||  secName == "pi-"        ||  secName == "pi0"        ||
             secName == "kaon+"  ||  secName == "kaon-"      ||  secName == "kaon0S"     ||
             secName == "kaon0L" ||  secName == "kaon0"      ||  secName == "anti_kaon0" ||
             secName == "eta"    ||  secName == "eta_prime"  ||  secName == "rho" ) {
          sumEnergySecondaries += secTE;
        } else if (  secName == "neutron"   ||  secName == "proton"  ||
                     secName == "deuteron"  ||  secName == "triton"  ||
                     secName == "alpha"     ||  secName == "He3" ) {	
          sumEnergySecondaries += secKE;
        } else if (  secName == "lambda"  ||  secName == "sigma0"  ||  secName == "sigma+"  ||
                     secName == "sigma-"  ||  secName == "xi-"     || secName == "xi0"      ||
                     secName == "omega-" ) {
          sumEnergySecondaries += secTE - 0.939;
        } else if ( secName == "anti_proton"  ||  secName == "anti_neutron" ||
                    secName == "anti_lambda"  ||  secName == "anti_sigma0"  ||
                    secName == "anti_sigma+"  ||  secName == "anti_sigma-"  ||
                    secName == "anti_xi-"     ||  secName == "anti_xi0"     ||
                    secName == "anti_omega-" ) {
          sumEnergySecondaries += secTE + 0.939;
        } else if ( secName == "electron" ) {
          sumEnergySecondaries += secKE;
        } else if ( secName == "positron" ) {
          sumEnergySecondaries += secTE + secMass;
        } else if ( secName == "mu-"     ||  secName == "mu+"          ||
                    secName == "tau-"    ||  secName == "tau+"         ||
                    secName == "nu_e"    ||  secName == "anti_nu_e"    ||
                    secName == "nu_mu"   ||  secName == "anti_nu_mu"   ||  
                    secName == "nu_tau"  ||  secName == "anti_nu_tau" ) { 
          sumEnergySecondaries += secTE;
        } else {  // It must be a fragment
          sumEnergySecondaries += secKE;
          //G4cout << " Other particle: " << secName << G4endl; //DEBUG
        }

        delete aChange->GetSecondary( j );

      }  // end loop on produced secondaries

      //G4cout << "\t event=" << loop << "\t sumEnergySecondaries=" 
      //       << sumEnergySecondaries << " GeV" << G4endl; //DEBUG

      const G4double thresholdOverProduction = 10.0;  // In GeV
      const G4double thresholdUnderProduction = 1.0;  // In GeV
      if ( sumEnergySecondaries - initialTotalEnergy > thresholdOverProduction ) {
        countOverProduction++;
        G4cerr << " ***VIOLATION-OVER*** event=" << loop 
               << " sumEnergySecondaries=" << sumEnergySecondaries << G4endl;
      } else if ( initialTotalEnergy - sumEnergySecondaries > thresholdUnderProduction ) {
        countUnderProduction++;
        G4cerr << " ***VIOLATION-UNDER*** event=" << loop 
               << " sumEnergySecondaries=" << sumEnergySecondaries << G4endl;
      }

      if ( sumEnergySecondaries > maxSumEnergySecondaries ) {
        maxSumEnergySecondaries = sumEnergySecondaries;
      } 
      if ( sumEnergySecondaries < minSumEnergySecondaries ) {
        minSumEnergySecondaries = sumEnergySecondaries;
      } 

      cumulatedSumEnergySecondaries += sumEnergySecondaries;
      cumulatedNumberSecondaries += nSec;

      aChange->Clear();
    }  // end event loop

    G4cout << "\t End event loop." << G4endl;

    G4double r = G4Random::getTheEngine()->flat();
    G4cout << "\t Final random number = " << r << G4endl;

    // Print some information at the end of the Run
    G4cout << "\t -------------  RUN SUMMARY  -------------" << G4endl
           << "\t average multiplicity = " 
           << (cumulatedNumberSecondaries*1.0)/(numberOfEvents*1.0) << G4endl
           << "\t countOverProduction  = " << countOverProduction << G4endl
           << "\t countUnderProduction = " << countUnderProduction << G4endl
           << "\t initialTotalEnergy           = " 
           << initialTotalEnergy << " GeV" << G4endl
           << "\t average sumEnergySecondaries = " 
           << cumulatedSumEnergySecondaries/(numberOfEvents*1.0) << " GeV" << G4endl
           << "\t maxSumEnergySecondaries      = " 
           << maxSumEnergySecondaries << " GeV" << G4endl
           << "\t minSumEnergySecondaries      = " 
           << minSumEnergySecondaries << " GeV" << G4endl
           << "\t -----------------------------------------" << G4endl
           << G4endl;

  } // End of the loop over the lines of the input file

  delete runManager;
  inputFile.close();
  return 0;
}
