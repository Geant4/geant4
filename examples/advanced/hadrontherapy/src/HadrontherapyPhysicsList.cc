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
// $Id: HadrontherapyPhysicsList.cc,v 1.0
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhysicsListMessenger.hh"
#include "HadrontherapyParticles.hh"
#include "Decay.hh"
#include "EMPhotonStandard.hh"
#include "EMPhotonEPDL.hh"
#include "EMPhotonPenelope.hh"
#include "EMElectronStandard.hh"
#include "EMElectronEEDL.hh"
#include "EMElectronPenelope.hh"
#include "EMPositronStandard.hh"
#include "EMPositronPenelope.hh"
#include "EMMuonStandard.hh"
#include "EMHadronIonLowEICRU49.hh"
#include "EMHadronIonLowEZiegler1977.hh"
#include "EMHadronIonLowEZiegler1985.hh"
#include "EMHadronIonStandard.hh"
#include "HIProtonNeutronPrecompound.hh"
#include "HIProtonNeutronBertini.hh"
#include "HIProtonNeutronBinary.hh"
#include "HIProtonNeutronLEP.hh"
#include "HIProtonNeutronPrecompoundGEM.hh"
#include "HIProtonNeutronPrecompoundFermi.hh"
#include "HIProtonNeutronPrecompoundGEMFermi.hh"
#include "HIPionBertini.hh"
#include "HIPionLEP.hh"
#include "HIIonLEP.hh"
#include "HEHadronIonLElastic.hh"
#include "HEHadronIonBertiniElastic.hh"
#include "HEHadronIonQElastic.hh"
#include "HEHadronIonUElastic.hh"
#include "HRMuonMinusCapture.hh"


HadrontherapyPhysicsList::HadrontherapyPhysicsList(): G4VModularPhysicsList(),
                                                      decayIsRegistered(false),
						      emElectronIsRegistered(false), 
						      emPositronIsRegistered(false), 
						      emPhotonIsRegistered(false), 
						      emIonIsRegistered(false),
						      emMuonIsRegistered(false),
						      hadrElasticHadronIonIsRegistered(false),
                                                      hadrInelasticPionIsRegistered(false),
                                                      hadrInelasticIonIsRegistered(false),
                                                      hadrInelasticProtonNeutronIsRegistered(false),
                                                      hadrAtRestMuonIsRegistered(false)
{
  // The secondary production threshold is set to 10. mm
  // for all the particles in all the experimental set-up
  // The phantom is defined as a Geant4 Region. Here the cut is fixed to 0.001 mm
  defaultCutValue = 0.01 * mm;

  // Messenger: it is possible to activate physics processes and models interactively 
  messenger = new HadrontherapyPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // Register all the particles involved in the experimental set-up
  RegisterPhysics( new HadrontherapyParticles("particles") );
}

HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{ 
  delete messenger;
}

void HadrontherapyPhysicsList::AddPhysicsList(const G4String& name)
{
  G4cout << "Adding PhysicsList component " << name << G4endl;
  

  // ****************
  // *** A. DECAY ***
  // ****************


  if (name == "Decay") 
    {
      if (decayIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new Decay(name) );
	  decayIsRegistered = true;
	}
    }


  // ***************************************
  // *** B. ELECTROMAGNETIC INTERACTIONS ***
  // ***************************************


  // ***************
  // *** Photons ***
  // ***************
  
  // *** Option I: Standard 

  if (name == "EM-Photon-Standard") 
    {
      if (emPhotonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMPhotonStandard(name) );
	  emPhotonIsRegistered = true;
	}
    }


  // *** Option II: Low Energy based on the Livermore libraries

  if (name == "EM-Photon-EPDL") 
    {
      if (emPhotonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing"
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMPhotonEPDL(name) );
	  emPhotonIsRegistered = true;
	}
    } 


  // *** Option III: Low Energy Penelope

  if (name == "EM-Photon-Penelope")
    {
      if (emPhotonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMPhotonPenelope(name) );
	  emPhotonIsRegistered = true;
	}
    }


  // *****************
  // *** Electrons ***
  // *****************
  
  // *** Option I: Standard

  if (name == "EM-Electron-Standard") 
    {
      if (emElectronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" 
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMElectronStandard(name) );	  
	  emElectronIsRegistered = true;
	}
    }


  // *** Option II: Low Energy based on the Livermore libraries

  if (name == "EM-Electron-EEDL") 
    {
      if (emElectronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing"                  
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMElectronEEDL(name) );
	  emElectronIsRegistered = true;
	}
    } 


  // *** Option III: Low Energy Penelope

  if (name == "EM-Electron-Penelope")
    {
      if (emElectronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- electron List already existing"                  
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMElectronPenelope(name) );
	  emElectronIsRegistered = true;
	}
    }


  // *****************
  // *** Positrons ***
  // *****************

  // *** Option I: Standard
  if (name == "EM-Positron-Standard") 
    {
      if (emPositronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing"                  
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMPositronStandard(name) );
	  emPositronIsRegistered = true;
	}
    }


  // *** Option II: Low Energy Penelope

  if (name == "EM-Positron-Penelope") 
    {
      if (emPositronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing"   
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMPositronPenelope(name) );
	  emPositronIsRegistered = true;
	}
    }
 

  // ************************
  // *** Hadrons and Ions ***
  // ************************

  // *** Option I: Low Energy with ICRU49 stopping power parametrisation
  
  if (name == "EM-HadronIon-LowE") 
    {
      if (emIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMHadronIonLowEICRU49(name) );
	  emIonIsRegistered = true;
	}
    }


  // *** Option II: Low Energy with Ziegler 1977 stopping power parametrisation

  if (name == "EM-HadronIon-LowEZiegler1977") 
    {
      if (emIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMHadronIonLowEZiegler1977(name) );
	  emIonIsRegistered = true;
	}
    }


  // *** Option III: Low Energy with Ziegler 1985 stopping power parametrisation

  if (name == "EM-HadronIon-LowEZiegler1985") 
    {
      if (emIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMHadronIonLowEZiegler1985(name) );
	  emIonIsRegistered = true;
	}
    }


  // *** Option IV: Standard

  if (name == "EM-HadronIon-Standard") 
    {
      if (emIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- ion List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMHadronIonStandard(name) );
	  emIonIsRegistered = true;
	}
    }


  // *************
  // *** Muons ***
  // *************

  // *** Option I: Standard 

  if (name == "EM-Muon-Standard") 
    {
      if (emMuonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new EMMuonStandard(name) );
	  emMuonIsRegistered = true;
	}
    }


  // ********************************
  // *** C. HADRONIC INTERACTIONS ***
  // ********************************


  // ******************************
  // *** C.1. ELASTIC PROCESSES ***
  // ******************************


  // ************************
  // *** Hadrons and Ions ***
  // ************************

  // *** Option I: GHEISHA like LEP model

  if (name == "HadronicEl-HadronIon-LElastic") 
    {
      if (hadrElasticHadronIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- hadronic Elastic Scattering List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HEHadronIonLElastic(name) );
	  hadrElasticHadronIonIsRegistered = true;
	}
    }
  

  // *** Option II: Bertini model

  if (name == "HadronicEl-HadronIon-Bert") 
    {
      if (hadrElasticHadronIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- hadronic Elastic Scattering List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HEHadronIonBertiniElastic(name) );
	  hadrElasticHadronIonIsRegistered = true;
	}
    }


  // *** Option III: Process G4QElastic

  if (name == "HadronicEl-HadronIon-QElastic") 
    {
      if (hadrElasticHadronIonIsRegistered) 
        {
          G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- hadronic Elastic Scattering List already existing" 
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HEHadronIonQElastic(name) );
          hadrElasticHadronIonIsRegistered = true;
	}
    }


  // *** Option III: Process G4UHadronElasticProcess

  if (name == "HadronicEl-HadronIon-UElastic") 
    {
     if (hadrElasticHadronIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- hadronic Elastic Scattering List already existing" 
		 << G4endl;
	} 
     else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
	         << " is registered" << G4endl;
	  RegisterPhysics( new HEHadronIonUElastic(name) );
	  hadrElasticHadronIonIsRegistered = true;
	}
    }
   
   
  // ********************************
  // *** C.2. INELASTIC PROCESSES ***
  // ********************************


  // *************
  // *** Pions ***
  // *************

  // *** Option I: Bertini model

  if (name == "HadronicInel-Pion-Bertini") 
    {
      if (hadrInelasticPionIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIPionBertini(name) );
	  hadrInelasticPionIsRegistered = true;
	}
    }


  // *** Option II: GHEISHA like LEP model

  if (name == "HadronicInel-Pion-LEP") 
    {
      if (hadrInelasticPionIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIPionLEP(name) );
	  hadrInelasticPionIsRegistered = true;
	}
    }


  // ************
  // *** Ions ***
  // ************

  // *** Option I: GHEISHA like LEP model

  if (name == "HadronicInel-Ion-LEP") 
    {
      if (hadrInelasticIonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIIonLEP(name) );
	  hadrInelasticIonIsRegistered = true;
	}
    }


  // *************************
  // *** Protons, Neutrons ***
  // *************************

  // *** Option I: GHEISHA like LEP model

  if (name == "HadronicInel-ProtonNeutron-LEP") 
    {
      if (hadrInelasticProtonNeutronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIProtonNeutronLEP(name) );
	  hadrInelasticProtonNeutronIsRegistered = true;
	}
    }


  // *** Option II: Bertini Cascade Model

  if (name == "HadronicInel-ProtonNeutron-Bert") 
    {
      if (hadrInelasticProtonNeutronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIProtonNeutronBertini(name) );
	  hadrInelasticProtonNeutronIsRegistered = true;
	}
    }


  // *** Option III: Binary Cascade Model

  if (name == "HadronicInel-ProtonNeutron-Bin") 
    {
      if (hadrInelasticProtonNeutronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIProtonNeutronBinary(name) );
	  hadrInelasticProtonNeutronIsRegistered = true;
	}
    }


  // *** Option IV: Precompound Model combined with Default Evaporation

  if (name == "HadronicInel-ProtonNeutron-Prec") 
    {
      if (hadrInelasticProtonNeutronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIProtonNeutronPrecompound(name) );
	  hadrInelasticProtonNeutronIsRegistered = true;
	}
    }


  // *** Option V: Precompound Model combined with GEM Evaporation

  if (name == "HadronicInel-ProtonNeutron-PrecGEM") 
    {
      if (hadrInelasticProtonNeutronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIProtonNeutronPrecompoundGEM(name) );
	  hadrInelasticProtonNeutronIsRegistered = true;
	}
    }


  // *** Option VI: Precompound Model combined with default Evaporation 
  //                and Fermi Break-up model

  if (name == "HadronicInel-ProtonNeutron-PrecFermi") 
    {
      if (hadrInelasticProtonNeutronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIProtonNeutronPrecompoundFermi(name) );
	  hadrInelasticProtonNeutronIsRegistered = true;
	}
    }


  // *** Option VII: Precompound Model combined with GEM Evaporation 
  //                 and Fermi Break-up model

  if (name == "HadronicInel-ProtonNeutron-PrecGEMFermi") 
    {
      if (hadrInelasticProtonNeutronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HIProtonNeutronPrecompoundGEMFermi(name) );
	  hadrInelasticProtonNeutronIsRegistered = true;
	}
    }   


  // ******************************
  // *** C.3. AT-REST PROCESSES ***
  // ******************************


  // **************
  // *** Muons- ***
  // **************

  // *** Option I: Muon Minus Capture at Rest

  if (name == "HadronicAtRest-MuonMinus-Capture") 
    {
      if (hadrAtRestMuonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HRMuonMinusCapture(name) );
	  hadrAtRestMuonIsRegistered = true;
	}
    }

}

void HadrontherapyPhysicsList::SetCuts()
{  
  // Set the threshold of production equal to the defaultCutValue
  // in the experimental set-up
  G4VUserPhysicsList::SetCutsWithDefault();
     
  G4double lowlimit=250*eV;
  G4ProductionCutsTable::GetProductionCutsTable() ->SetEnergyRange(lowlimit, 100.*GeV);
  // Definition of a smaller threshold of production in the phantom region
  // where high accuracy is required in the energy deposit calculation

  G4String regionName = "PhantomLog";
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);
  G4ProductionCuts* cuts = new G4ProductionCuts ;
  G4double regionCut = 0.01*mm;
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("gamma"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e-"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e+"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("proton"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("genericIons"));
  region -> SetProductionCuts(cuts);

  if (verboseLevel>0) DumpCutValuesTable();
}


