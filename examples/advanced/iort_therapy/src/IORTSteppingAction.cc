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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "IORTSteppingAction.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UserEventAction.hh"  // NOT INCLUDED IN geant4.9.3p01_version

#include "IORTAnalysisManager.hh"

#include "IORTRunAction.hh"

/////////////////////////////////////////////////////////////////////////////
IORTSteppingAction::IORTSteppingAction( IORTRunAction *run)
{
    runAction = run;
}

/////////////////////////////////////////////////////////////////////////////
IORTSteppingAction::~IORTSteppingAction()
{
}

/////////////////////////////////////////////////////////////////////////////
void IORTSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
    /*
    // USEFULL METHODS TO RETRIEVE INFORMATION DURING THE STEPS
    if( (aStep->GetTrack()->GetVolume()->GetName() == "DetectorPhys") 
    && aStep->GetTrack()->GetDefinition()->GetParticleName() == "proton")
    //G4int evtNb = G4RunManager::GetRunManager()->GetCurrentEvent() -> GetEventID();
    {
    G4cout << "ENERGIA:   " << aStep->GetTrack()->GetKineticEnergy() 
    << "   VOLUME  " << aStep->GetTrack()->GetVolume()->GetName()
    << "   MATERIALE    " <<  aStep -> GetTrack() -> GetMaterial() -> GetName()
    << "   EVENTO     " << G4RunManager::GetRunManager()->GetCurrentEvent() -> GetEventID()
    << "   POS     " << aStep->GetTrack()->GetPosition().x()
    << G4endl;
    }
    */

    if( aStep->GetTrack()->GetVolume()->GetName() == "NewDetectorPhys"){
#ifdef G4ANALYSIS_USE_ROOT
	G4ParticleDefinition *def = aStep->GetTrack()->GetDefinition();
	G4double secondaryParticleKineticEnergy =  aStep->GetTrack()->GetKineticEnergy();     
	G4String particleType = def->GetParticleType(); // particle type = nucleus for d, t, He3, alpha, and heavier nuclei
	G4String particleName = def->GetParticleName(); // e.g. for alpha: the name = "alpha" and type = "nucleus"
	if(particleType == "nucleus") {
	    G4int A = def->GetBaryonNumber();
	    G4double Z = def->GetPDGCharge();
	    G4double posX = aStep->GetTrack()->GetPosition().x() / CLHEP::cm;
	    G4double posY = aStep->GetTrack()->GetPosition().y() / CLHEP::cm;
	    G4double posZ = aStep->GetTrack()->GetPosition().z() / CLHEP::cm;
	    G4double energy = secondaryParticleKineticEnergy / A / CLHEP::MeV;

	    IORTAnalysisManager* analysisMgr =  IORTAnalysisManager::GetInstance();   
	    analysisMgr->FillFragmentTuple(A, Z, energy, posX, posY, posZ);
	} else if(particleName == "proton") {   // proton (hydrogen-1) is a special case
	    G4double posX = aStep->GetTrack()->GetPosition().x() / CLHEP::cm ;
	    G4double posY = aStep->GetTrack()->GetPosition().y() / CLHEP::cm ;
	    G4double posZ = aStep->GetTrack()->GetPosition().z() / CLHEP::cm ;
	    G4double energy = secondaryParticleKineticEnergy * CLHEP::MeV;    // Hydrogen-1: A = 1, Z = 1
	    IORTAnalysisManager::GetInstance()->FillFragmentTuple(1, 1.0, energy, posX, posY, posZ);
	}

	G4String secondaryParticleName =  def -> GetParticleName();  
	//G4cout <<"Particle: " << secondaryParticleName << G4endl;
	//G4cout <<"Energy: " << secondaryParticleKineticEnergy << G4endl;
	IORTAnalysisManager* analysis =  IORTAnalysisManager::GetInstance();   
	//There is a bunch of stuff recorded with the energy 0, something should perhaps be done about this.
	if(secondaryParticleName == "proton") {
	    analysis->hydrogenEnergy(secondaryParticleKineticEnergy / CLHEP::MeV);
	}
	if(secondaryParticleName == "deuteron") {
	    analysis->hydrogenEnergy((secondaryParticleKineticEnergy/2) / CLHEP::MeV);
	}
	if(secondaryParticleName == "triton") {
	    analysis->hydrogenEnergy((secondaryParticleKineticEnergy/3) / CLHEP::MeV);
	}
	if(secondaryParticleName == "alpha") {
	    analysis->heliumEnergy((secondaryParticleKineticEnergy/4) / CLHEP::MeV);
	}
	if(secondaryParticleName == "He3"){
	    analysis->heliumEnergy((secondaryParticleKineticEnergy/3) / CLHEP::MeV);		
	}
#endif

	aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    }

    // Electromagnetic and hadronic processes of primary particles in the phantom
    //setting phantomPhys correctly will break something here fixme
    if ((aStep -> GetTrack() -> GetTrackID() == 1) &&
	    (aStep -> GetTrack() -> GetVolume() -> GetName() == "PhantomPhys") &&
	    (aStep -> GetPostStepPoint() -> GetProcessDefinedStep() != NULL))
    {
	G4String process = aStep -> GetPostStepPoint() -> 
	    GetProcessDefinedStep() -> GetProcessName();

	if ((process == "Transportation") || (process == "StepLimiter")) {;}
	else 
	{
	    if ((process == "msc") || (process == "hLowEIoni") || (process == "hIoni")) 
	    { 
		runAction -> AddEMProcess();
	    } 
	    else  
	    {
		runAction -> AddHadronicProcess();

		if ( (process != "LElastic") && (process != "ProtonInelastic") && (process != "hElastic") )
		    G4cout << "Warning! Unknown proton process: "<< process << G4endl;
	    }
	}         
    }

    // Retrieve information about the secondary particles originated in the phantom

    G4SteppingManager*  steppingManager = fpSteppingManager;

    // check if it is alive
    //if(theTrack-> GetTrackStatus() == fAlive) { return; }

    // Retrieve the secondary particles
    G4TrackVector* fSecondary = steppingManager -> GetfSecondary();

    for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
    { 
	G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); 

	if (volumeName == "phantomPhys")
	{
#ifdef G4ANALYSIS_USE_ROOT   
	    G4String secondaryParticleName =  (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();  
	    G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy();     

	    IORTAnalysisManager* analysis =  IORTAnalysisManager::GetInstance();   

	    if (secondaryParticleName == "e-")
		analysis -> electronEnergyDistribution(secondaryParticleKineticEnergy/CLHEP::MeV);

	    if (secondaryParticleName == "gamma")
		analysis -> gammaEnergyDistribution(secondaryParticleKineticEnergy/CLHEP::MeV);

	    if (secondaryParticleName == "deuteron")
		analysis -> deuteronEnergyDistribution(secondaryParticleKineticEnergy/CLHEP::MeV);

	    if (secondaryParticleName == "triton")
		analysis -> tritonEnergyDistribution(secondaryParticleKineticEnergy/CLHEP::MeV);

	    if (secondaryParticleName == "alpha")
		analysis -> alphaEnergyDistribution(secondaryParticleKineticEnergy/CLHEP::MeV);

	    G4double z = (*fSecondary)[lp1]-> GetDynamicParticle() -> GetDefinition() -> GetPDGCharge();
	    if (z > 0.)
	    {      
		G4int a = (*fSecondary)[lp1]-> GetDynamicParticle() -> GetDefinition() -> GetBaryonNumber();
		G4int electronOccupancy = (*fSecondary)[lp1] ->  GetDynamicParticle() -> GetTotalOccupancy(); 

		// If a generic ion is originated in the detector, its baryonic number, PDG charge, 
		// total number of electrons in the orbitals are stored in a ntuple 
		analysis -> genericIonInformation(a, z, electronOccupancy, secondaryParticleKineticEnergy/CLHEP::MeV);
	    }
#endif
	}
    }
}
