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
// $Id: RemSimSteppingAction.cc,v 1.10 2006-06-29 16:24:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "RemSimSteppingAction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "G4TrackStatus.hh"
#include "RemSimSteppingActionMessenger.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

RemSimSteppingAction::RemSimSteppingAction(RemSimPrimaryGeneratorAction* primary):
  primaryAction(primary)
{ 
  hadronic = "Off";
  messenger = new RemSimSteppingActionMessenger(this);
}

RemSimSteppingAction::~RemSimSteppingAction()
{ 
  delete messenger;
}

void RemSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4String oldVolumeName = aStep -> GetPreStepPoint()->  
                                      GetPhysicalVolume()-> GetName();  
 
  // Store in histograms useful information concerning primary particles

#ifdef G4ANALYSIS_USE    
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();

  G4VPhysicalVolume* volume = aStep -> GetPostStepPoint() -> GetPhysicalVolume();
 
 
  // Particle reaching the astronaut
  if(oldVolumeName != "phantom") 
    { 
      if (volume) 
	{
	  G4String volumeName = volume -> GetName();
	 G4String particleName = (aStep -> GetTrack() -> GetDynamicParticle()
			   -> GetDefinition() -> GetParticleName());

	 if (volumeName == "phantom") 
	    { 
	      G4double particleEnergy = aStep -> GetTrack() -> GetKineticEnergy();
	     
	    if(aStep -> GetTrack() -> GetTrackID()== 1) // primary particles
		{ 
	         G4double initialEnergy = primaryAction -> GetInitialEnergy(); 
	        
	      if ( particleName == "proton" || 
                   particleName == "alpha"  ||
		   particleName == "IonO16" || 
		   particleName == "IonC12" || 
                   particleName == "IonFe52"|| 
		   particleName == "IonSi28")
		{
		 G4int baryon = aStep -> GetTrack() -> GetDefinition() -> GetBaryonNumber();
		
                 // Initial energy of primary particles impinging on the phantom
		 analysis -> PrimaryInitialEnergyIn((initialEnergy/baryon)/MeV);
                 
		 // Energy of primary particles impinging on the phantom
	         analysis -> PrimaryEnergyIn((particleEnergy/baryon)/MeV);
		}
		}
	        // secondary particle reaching the astronaut
	        else {
       
		  // if i =0  secondary proton, i =1 neutron, i=2 pion, i=3 alpha, 
		  // i =4 other, i=5 electron, i = 6 gamma, i=7 positrons, 
		  // i=8 muons, i=9 neutrinos

		if (particleName == "proton") 
		  {
		    analysis -> SecondaryReachingThePhantom(0);
		    analysis -> SecondaryProtonReachingThePhantom(particleEnergy/MeV); 
		  }
	  else
	    {
	      if (particleName == "neutron") 
		{
		  analysis -> SecondaryReachingThePhantom(1);
	          analysis -> SecondaryNeutronReachingThePhantom(particleEnergy/MeV);
		}
		  else{
                     if (particleName == "pi+" ||
			 particleName == "pi-" ||particleName == "pi0" ) 
		       {
			 analysis -> SecondaryReachingThePhantom(2);
			 analysis -> SecondaryPionReachingThePhantom(particleEnergy/MeV);
		       }
		     else{
		       if (particleName == "alpha") 
			 {
			   analysis -> SecondaryReachingThePhantom(3);
			   analysis -> SecondaryAlphaReachingThePhantom(particleEnergy/MeV);
			 }
		       else{
			 if (particleName == "e+") 
			   {analysis -> SecondaryReachingThePhantom(7);
			   analysis -> SecondaryPositronReachingThePhantom(particleEnergy/MeV);
			   } 
			 else{
			   if (particleName == "e-") 
			     {analysis -> SecondaryReachingThePhantom(5);
			     analysis -> SecondaryElectronReachingThePhantom(particleEnergy/MeV);
			     }			    
			   else{
			     if (particleName == "gamma") 
			       {analysis -> SecondaryReachingThePhantom(6); 
			       analysis -> SecondaryGammaReachingThePhantom(particleEnergy/MeV);	
			       }	
			     else{
			       if (particleName == "mu+" 
				   ||particleName == "mu-" ) 
				 {analysis -> SecondaryReachingThePhantom(8);
				 analysis -> SecondaryMuonReachingThePhantom(particleEnergy/MeV);
				 }
			       else {
				 if (particleName == "nu_e" || particleName == "nu_mu" ||
                                     particleName == "anti_nu_e" || particleName == "anti_nu_mu")
				 { analysis -> SecondaryReachingThePhantom(9);
				}
			       else{ 
				 analysis -> SecondaryReachingThePhantom(4);
				 analysis -> SecondaryOtherReachingThePhantom(particleEnergy/MeV);
			       }
			       }
			     }
			   }
			 }
		       }
		     }
		  }
	         }
	      }
	    }
	}
    }
   
	   
  // Primar particles outgoing the phantom
  if(oldVolumeName == "phantom") 
    { 
      if (volume) 
	{
	  G4String volumeName = volume -> GetName();
          G4String particleName = (aStep -> GetTrack() -> GetDynamicParticle()
			   -> GetDefinition() -> GetParticleName()); 
	  if (volumeName != "phantom") 
	    {
	      if(aStep -> GetTrack() -> GetTrackID()== 1) // primary particles
		{
		  G4double initialEnergy = primaryAction -> GetInitialEnergy(); 
		  G4double particleEnergy = aStep->GetTrack() -> GetKineticEnergy();
		  if ( particleName == "proton" || 
		       particleName == "alpha"  || 
		       particleName == "IonO16" || 
		       particleName == "IonC12" || 
		       particleName == "IonFe52"|| 
		       particleName == "IonSi28" )
		    { 
		      G4int baryon = aStep -> GetTrack() -> GetDefinition() -> GetBaryonNumber();
		      // plot of MeV/nucl
		      analysis -> PrimaryInitialEnergyOut((initialEnergy/baryon)/MeV); 
		      analysis -> PrimaryEnergyOut((particleEnergy/baryon)/MeV);
		    }
		}
	    }
	}
    }

  // Analysis of the secondary particles generated in the phantom

  G4SteppingManager*  steppingManager = fpSteppingManager;
  G4Track* theTrack = aStep -> GetTrack();

  // check if it is alive
  if(theTrack-> GetTrackStatus() == fAlive) { return; }

  // Retrieve the secondary particles
   G4TrackVector* fSecondary = steppingManager -> GetfSecondary();
     
   for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
     { 
       // Retrieve the info about the generation of secondary particles in the phantom and
       // in the vehicle
       G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); 
       G4String secondaryParticleName =  (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();  
       G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy(); 
       G4String process = (*fSecondary)[lp1]-> GetCreatorProcess()-> GetProcessName();   
      
       // If the secondaries are originated in the phantom....
       if (volumeName == "phantom")
	 {
	 
	   G4double slice = (*fSecondary)[lp1] -> GetPosition().z() ;
	   //G4cout << "Secondary particle in phantom: " << secondaryParticleName
	   //	  << "energy (MeV): " << secondaryParticleKineticEnergy << G4endl; 
	   //      << "creation slice: " << slice << G4endl; 
	   // if i =0  secondary proton, i =1 neutron, i=2 pion, i=3 alpha, 
	   // i =4 other, i=5 electron, i = 6 gamma, i=7 positrons, 
	   // i=8 muons, i=9 neutrinos
	   G4double translation = 15. * cm;
	   if (secondaryParticleName == "proton") 
	     {
	       analysis -> SecondaryInPhantom(0);
	       analysis -> SecondaryProtonInPhantom(secondaryParticleKineticEnergy/MeV); 
	       slice = slice + translation;
	       analysis -> SecondaryProtonInPhantomSlice(slice/cm); }
	   else
	     {
	       if (secondaryParticleName == "neutron") 
		 {
		   analysis -> SecondaryInPhantom(1);
		   analysis -> SecondaryNeutronInPhantom(secondaryParticleKineticEnergy/MeV);
		   slice = slice + translation;
		   analysis -> SecondaryNeutronInPhantomSlice(slice/cm); 
		 }
	       else{
		 if (secondaryParticleName == "pi+" ||
		     secondaryParticleName == "pi-"||secondaryParticleName == "pi0") 
		   {
		     analysis -> SecondaryInPhantom(2);
		     analysis -> SecondaryPionInPhantom(secondaryParticleKineticEnergy/MeV);
		     slice = slice + translation;
		     analysis -> SecondaryPionInPhantomSlice(slice/cm);  
		   }
		 else{
		   if (secondaryParticleName == "alpha") 
		     {
		       analysis -> SecondaryInPhantom(3);
		       analysis -> SecondaryAlphaInPhantom(secondaryParticleKineticEnergy/MeV);
		       slice = slice + translation;
		       analysis -> SecondaryAlphaInPhantomSlice(slice/cm);
		     }
		   else{
		     if (secondaryParticleName == "e+") 
		       {analysis -> SecondaryInPhantom(7);
		       analysis -> SecondaryPositronInPhantom(secondaryParticleKineticEnergy/MeV);
		       slice = slice + translation;
		       analysis -> SecondaryPositronInPhantomSlice(slice/cm);
		       } 
		     else{
		       if (secondaryParticleName == "e-") 
			 {analysis -> SecondaryInPhantom(5);
			 analysis -> SecondaryElectronInPhantom(secondaryParticleKineticEnergy/MeV);
			 slice = slice + translation;
			 analysis -> SecondaryElectronInPhantomSlice(slice/cm);
			 }			    
		       else{
			 if (secondaryParticleName == "gamma") 
			   {analysis -> SecondaryInPhantom(6); 
			   analysis -> SecondaryGammaInPhantom(secondaryParticleKineticEnergy/MeV);
			   slice = slice + translation;
			   analysis -> SecondaryGammaInPhantomSlice(slice/cm);
			   }	
			 else{
			   if (secondaryParticleName == "mu+" 
			       ||secondaryParticleName == "mu-" ) 
			     {analysis -> SecondaryInPhantom(8);
			     analysis -> SecondaryMuonInPhantom(secondaryParticleKineticEnergy/MeV);
			     slice = slice + translation;
			     analysis -> SecondaryMuonInPhantomSlice(slice/cm);
			     }
			   else{ 
			     if (secondaryParticleName == "nu_e" || secondaryParticleName == "nu_mu"
				 || secondaryParticleName == "anti_nu_e" || secondaryParticleName == "anti_nu_mu")
			       {  
				 analysis -> SecondaryInPhantom(9);
				 analysis -> SecondaryNeutrinoInPhantom(secondaryParticleKineticEnergy/MeV);} 
			     else{				  
			       analysis -> SecondaryInPhantom(4);
			       analysis -> SecondaryOtherInPhantom(secondaryParticleKineticEnergy/MeV);
			       slice = slice + translation;
			       analysis -> SecondaryOtherInPhantomSlice(slice/cm);
			     }
			   }
			 }
		       }
		     }
		   }

		 }
	       }
	     }
	 }

       // secondary particles produced in the multilayer + shielding
       else{ 
	 if (volumeName != "world")
	   {
	     if (secondaryParticleName == "proton") analysis -> SecondaryInVehicle(0);
	     else
	       {
		 if (secondaryParticleName == "neutron")  analysis -> SecondaryInVehicle(1);
		 else{
		   if (secondaryParticleName == "pi+" ||
		       secondaryParticleName == "pi-"||secondaryParticleName == "pi0")  analysis -> SecondaryInVehicle(2);
		   else{
		     if (secondaryParticleName == "alpha") analysis -> SecondaryInVehicle(3);
                         			 
		     else{
		       if (secondaryParticleName == "e+") analysis -> SecondaryInVehicle(7);
			     
		       else{
			 if (secondaryParticleName == "e-") analysis -> SecondaryInVehicle(5);
			 else{
			   if (secondaryParticleName == "gamma") analysis -> SecondaryInVehicle(6); 
			   else{
			     if (secondaryParticleName == "mu+" 
				 ||secondaryParticleName == "mu-" ) analysis -> SecondaryInVehicle(8);
			     else{ 
			       if (secondaryParticleName == "nu_e"|| secondaryParticleName == "nu_mu" ||
				   secondaryParticleName == "anti_nu_e" || secondaryParticleName == "anti_nu_mu")
				 {analysis -> SecondaryInVehicle(9);}
				
			       else{				  
				 analysis -> SecondaryInVehicle(4);
			       }
			     }
			   }
			 }
		       }
		     }
		   }
		 }
	       }
	   }
       }
     }
    
#endif

   //analysis of hadronic physics
   if (hadronic == "On")
     {
       if(aStep -> GetPostStepPoint() -> GetProcessDefinedStep() != NULL)
	 {
	   G4String process = aStep -> GetPostStepPoint() ->
	     GetProcessDefinedStep() ->GetProcessName();

	   if ((process != "Transportation") &&
	       (process != "msc") && 
	       (process != "LowEnergyIoni") &&
	       (process != "LowEnergyBrem") && 
	       (process != "eIoni") &&
	       (process != "hIoni") &&
	       (process != "eBrem") && 
	       (process != "compt") &&
	       (process != "phot")  && 
	       (process != "conv")  &&
	       (process != "annihil") &&
	       (process != "hLowEIoni") &&
	       (process != "LowEnBrem") && 
	       (process != "LowEnCompton") && 
	       (process != "LowEnPhotoElec") && 
	       (process != "LowEnRayleigh") && 
	       (process != "LowEnConversion"))
     	     G4cout << "Hadronic Process:" << process << G4endl;
	 }
     }
}
void RemSimSteppingAction::SetHadronicAnalysis(G4String value)    
{
  hadronic = value;
}
