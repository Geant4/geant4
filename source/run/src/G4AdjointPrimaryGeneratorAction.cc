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
// $Id: G4AdjointPrimaryGeneratorAction.cc 102435 2017-01-27 08:28:15Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointPrimaryGeneratorAction
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////

#include "G4AdjointPrimaryGeneratorAction.hh"

#include "G4PhysicalConstants.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh" 
#include "G4AdjointSimManager.hh" 
#include "G4AdjointPrimaryGenerator.hh"
#include "G4Gamma.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPrimaryGeneratorAction::G4AdjointPrimaryGeneratorAction()
  : Emin(0.), Emax(0.), EminIon(0.), EmaxIon(0.),
    index_particle(100000),
    radius_spherical_source(0.), fwd_ion(0), adj_ion(0), 
    ion_name("not_defined")
{
  theAdjointPrimaryGenerator= new G4AdjointPrimaryGenerator();

  PrimariesConsideredInAdjointSim[G4String("e-")]=false;
  PrimariesConsideredInAdjointSim[G4String("gamma")]=false;
  PrimariesConsideredInAdjointSim[G4String("proton")]=false;
  PrimariesConsideredInAdjointSim[G4String("ion")]=false;

  ListOfPrimaryFwdParticles.clear();
  ListOfPrimaryAdjParticles.clear();
  nb_fwd_gammas_per_event = 1;
  nb_adj_primary_gammas_per_event = 1;
  nb_adj_primary_electrons_per_event = 1;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPrimaryGeneratorAction::~G4AdjointPrimaryGeneratorAction()
{
  delete theAdjointPrimaryGenerator;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

   G4int evt_id=anEvent->GetEventID();
   size_t n=ListOfPrimaryAdjParticles.size();
   index_particle=size_t(evt_id)-n*(size_t(evt_id)/n);


   G4double E1=Emin;
   G4double E2=Emax;
   if (!ListOfPrimaryAdjParticles[index_particle]) UpdateListOfPrimaryParticles();//ion has not been created yet

   if (ListOfPrimaryAdjParticles[index_particle]->GetParticleName() == "adj_proton") {
     E1=EminIon;
     E2=EmaxIon;
   }
   if (ListOfPrimaryAdjParticles[index_particle]->GetParticleType() == "adjoint_nucleus") {
     G4int A= ListOfPrimaryAdjParticles[index_particle]->GetAtomicMass();
     E1=EminIon*A;
     E2=EmaxIon*A;
   }
   //Generate first the forwrad primaries
   theAdjointPrimaryGenerator->GenerateFwdPrimaryVertex(anEvent,ListOfPrimaryFwdParticles[index_particle],E1,E2);
   G4PrimaryVertex* fwdPrimVertex = anEvent->GetPrimaryVertex();

   p=fwdPrimVertex->GetPrimary()->GetMomentum();
   pos=fwdPrimVertex->GetPosition();
   G4double pmag=p.mag();
   G4double m0=ListOfPrimaryFwdParticles[index_particle]->GetPDGMass();
   G4double ekin=std::sqrt( m0*m0 + pmag*pmag) -m0;

   G4double weight_correction=1.;
   //For gamma generate the particle along the backward ray
   G4ThreeVector dir=-p/p.mag();

   /*if (ListOfPrimaryAdjParticles[index_particle]->GetParticleName() == "adj_gamma"){

      theAdjointPrimaryGenerator
              ->ComputeAccumulatedDepthVectorAlongBackRay(pos,dir,ekin,ListOfPrimaryAdjParticles[index_particle]);


      G4double distance = theAdjointPrimaryGenerator
       ->SampleDistanceAlongBackRayAndComputeWeightCorrection(weight_correction);

      //pos=pos+dir*distance;
      //fwdPrimVertex->SetPosition(pos[0],pos[1],pos[2]);
   }
   */
   weight_correction=1.;

   if (ListOfPrimaryFwdParticles[index_particle]  ==G4Gamma::Gamma() && nb_fwd_gammas_per_event>1 ){
	   G4double weight = (1./nb_fwd_gammas_per_event);
	   fwdPrimVertex->SetWeight(weight);
	   for (int i=0;i<nb_fwd_gammas_per_event-1;i++){
         G4PrimaryVertex* newFwdPrimVertex = new G4PrimaryVertex();
         newFwdPrimVertex->SetPosition(pos.x(),pos.y(),pos.z());
         newFwdPrimVertex->SetT0(0.);
         G4PrimaryParticle* aPrimParticle = new G4PrimaryParticle(ListOfPrimaryFwdParticles[index_particle],
               p.x(),p.y(),p.z());
         newFwdPrimVertex->SetPrimary(aPrimParticle);
         newFwdPrimVertex->SetWeight(weight);
         anEvent->AddPrimaryVertex(newFwdPrimVertex);
       }
   }

   //Now generate the adjoint primaries
   G4PrimaryVertex* adjPrimVertex = new G4PrimaryVertex();
   adjPrimVertex->SetPosition(pos.x(),pos.y(),pos.z());
   adjPrimVertex->SetT0(0.);
   G4PrimaryParticle* aPrimParticle = new G4PrimaryParticle(ListOfPrimaryAdjParticles[index_particle],
                         -p.x(),-p.y(),-p.z());

   adjPrimVertex->SetPrimary(aPrimParticle);
   anEvent->AddPrimaryVertex(adjPrimVertex);

   //The factor pi is to normalise the weight to the directional flux
   G4double adjoint_source_area = G4AdjointSimManager::GetInstance()->GetAdjointSourceArea();
   G4double adjoint_weight = weight_correction*ComputeEnergyDistWeight(ekin,E1,E2)*adjoint_source_area*pi;
   //if (ListOfPrimaryFwdParticles[index_particle]  ==G4Gamma::Gamma()) adjoint_weight = adjoint_weight/3.;
   if (ListOfPrimaryAdjParticles[index_particle]->GetParticleName() == "adj_gamma") {
	   //The weight will be corrected at the end of the track if splitted tracks are used
	   adjoint_weight = adjoint_weight/nb_adj_primary_gammas_per_event;
	   for (int i=0;i<nb_adj_primary_gammas_per_event-1;i++){
          G4PrimaryVertex* newAdjPrimVertex = new G4PrimaryVertex();
	      newAdjPrimVertex->SetPosition(pos.x(),pos.y(),pos.z());
	      newAdjPrimVertex->SetT0(0.);
	      aPrimParticle = new G4PrimaryParticle(ListOfPrimaryAdjParticles[index_particle],
	                              -p.x(),-p.y(),-p.z());
	      newAdjPrimVertex->SetPrimary(aPrimParticle);
	      newAdjPrimVertex->SetWeight(adjoint_weight);
	      anEvent->AddPrimaryVertex(newAdjPrimVertex);
	   }
   }
   else if  (ListOfPrimaryAdjParticles[index_particle]->GetParticleName() == "adj_electron") {
	   //The weight will be corrected at the end of the track if splitted tracks are used
	   adjoint_weight = adjoint_weight/nb_adj_primary_electrons_per_event;
	   for (int i=0;i<nb_adj_primary_electrons_per_event-1;i++){
          G4PrimaryVertex* newAdjPrimVertex = new G4PrimaryVertex();
	      newAdjPrimVertex->SetPosition(pos.x(),pos.y(),pos.z());
	      newAdjPrimVertex->SetT0(0.);
	      aPrimParticle = new G4PrimaryParticle(ListOfPrimaryAdjParticles[index_particle],
	                              -p.x(),-p.y(),-p.z());
	      newAdjPrimVertex->SetPrimary(aPrimParticle);
	      newAdjPrimVertex->SetWeight(adjoint_weight);
	      anEvent->AddPrimaryVertex(newAdjPrimVertex);
	   }
   }
   adjPrimVertex->SetWeight(adjoint_weight);

   //Call some methods of G4AdjointSimManager
   G4AdjointSimManager::GetInstance()->SetAdjointTrackingMode(true);
   G4AdjointSimManager::GetInstance()->ClearEndOfAdjointTrackInfoVectors();
   G4AdjointSimManager::GetInstance()->ResetDidOneAdjPartReachExtSourceDuringEvent();

  /* if ( !last_generated_part_was_adjoint ) {
	 
	 index_particle++;
	 if (index_particle >= ListOfPrimaryAdjParticles.size()) index_particle =0;
	 
	
	 G4double E1=Emin;
	 G4double E2=Emax;
	 if (!ListOfPrimaryAdjParticles[index_particle]) UpdateListOfPrimaryParticles();//ion has not been created yet
	
	 if (ListOfPrimaryAdjParticles[index_particle]->GetParticleName() == "adj_proton") {
		E1=EminIon;
		E2=EmaxIon;
	 }
	 if (ListOfPrimaryAdjParticles[index_particle]->GetParticleType() == "adjoint_nucleus") {
		G4int A= ListOfPrimaryAdjParticles[index_particle]->GetAtomicMass();
		E1=EminIon*A;
		E2=EmaxIon*A;
	 }
	 theAdjointPrimaryGenerator->GenerateAdjointPrimaryVertex(anEvent,
	 							  ListOfPrimaryAdjParticles[index_particle],
								  E1,E2);
	 G4PrimaryVertex* aPrimVertex = anEvent->GetPrimaryVertex();
  
  
  	  p=aPrimVertex->GetPrimary()->GetMomentum();
	  pos=aPrimVertex->GetPosition();
	  G4double pmag=p.mag();
	  
	  G4double m0=ListOfPrimaryAdjParticles[index_particle]->GetPDGMass();
	  G4double ekin=std::sqrt( m0*m0 + pmag*pmag) -m0;
  	
  
  	  //The factor pi is to normalise the weight to the directional flux
	  G4double adjoint_source_area = G4AdjointSimManager::GetInstance()->GetAdjointSourceArea();
  	  G4double adjoint_weight = ComputeEnergyDistWeight(ekin,E1,E2)*adjoint_source_area*pi;

  	  aPrimVertex->SetWeight(adjoint_weight);

	  last_generated_part_was_adjoint =true;
	  G4AdjointSimManager::GetInstance()->SetAdjointTrackingMode(true);
	  G4AdjointSimManager::GetInstance()->RegisterAdjointPrimaryWeight(adjoint_weight);
   }
   else {
          //fwd particle equivalent to the last generated adjoint particle ios generated	
  	  G4PrimaryVertex* aPrimVertex = new G4PrimaryVertex();
	  aPrimVertex->SetPosition(pos.x(),pos.y(),pos.z());
  	  aPrimVertex->SetT0(0.);
	  G4PrimaryParticle* aPrimParticle = new G4PrimaryParticle(ListOfPrimaryFwdParticles[index_particle],
  							   -p.x(),-p.y(),-p.z()); 
	
	  aPrimVertex->SetPrimary(aPrimParticle);
 	  anEvent->AddPrimaryVertex(aPrimVertex); 						   
	  last_generated_part_was_adjoint =false;
	  G4AdjointSimManager::GetInstance()->SetAdjointTrackingMode(false);
	*/

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEmin(G4double val)
{
  Emin=val;
  EminIon=val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEmax(G4double val)
{
  Emax=val;
  EmaxIon=val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEminIon(G4double val)
{
  EminIon=val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEmaxIon(G4double val)
{
  EmaxIon=val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointPrimaryGeneratorAction::ComputeEnergyDistWeight(G4double E ,G4double E1, G4double E2)
{
 //  We generate N numbers of primaries with  a 1/E energy law distribution.
 //  We have therefore  an energy distribution function    
 // 		f(E)=C/E  (1)
 //  with   C a constant that is such that
 //		N=Integral(f(E),E1,E2)=C.std::log(E2/E1)  (2)
 //  Therefore from (2) we get		
 // 		C=N/ std::log(E2/E1) (3) 
 //  and        
 //		f(E)=N/ std::log(E2/E1)/E (4)
 //For the adjoint simulation we need a energy distribution  f'(E)=1..
 //To get that we need therefore to apply a weight to the primary
 //		W=1/f(E)=E*std::log(E2/E1)/N  
 //
  return std::log(E2/E1)*E/G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetSphericalAdjointPrimarySource(G4double radius, G4ThreeVector center_pos)
{ 
  radius_spherical_source = radius;
  center_spherical_source = center_pos;
  type_of_adjoint_source ="Spherical";
  theAdjointPrimaryGenerator->SetSphericalAdjointPrimarySource(radius,center_pos);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetAdjointPrimarySourceOnAnExtSurfaceOfAVolume(const G4String& volume_name)
{
  type_of_adjoint_source ="ExternalSurfaceOfAVolume";
  theAdjointPrimaryGenerator->SetAdjointPrimarySourceOnAnExtSurfaceOfAVolume(volume_name);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::ConsiderParticleAsPrimary(const G4String& particle_name)
{
  if (PrimariesConsideredInAdjointSim.find(particle_name) != PrimariesConsideredInAdjointSim.end()){
  	PrimariesConsideredInAdjointSim[particle_name]=true;
  }
  UpdateListOfPrimaryParticles();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::NeglectParticleAsPrimary(const G4String& particle_name)
{
  if (PrimariesConsideredInAdjointSim.find(particle_name) != PrimariesConsideredInAdjointSim.end()){
  	PrimariesConsideredInAdjointSim[particle_name]= false;
  }
  UpdateListOfPrimaryParticles();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::UpdateListOfPrimaryParticles()
{
    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    ListOfPrimaryFwdParticles.clear();
    ListOfPrimaryAdjParticles.clear();
    std::map<G4String, G4bool>::iterator iter;
    for( iter = PrimariesConsideredInAdjointSim.begin(); iter != PrimariesConsideredInAdjointSim.end(); ++iter ) {
 	if(iter->second) {
		G4String fwd_particle_name = iter->first;
		if ( fwd_particle_name != "ion") {
			G4String adj_particle_name = G4String("adj_") + fwd_particle_name;
			ListOfPrimaryFwdParticles.push_back(theParticleTable->FindParticle(fwd_particle_name));
			ListOfPrimaryAdjParticles.push_back(theParticleTable->FindParticle(adj_particle_name));
			/*
			if ( fwd_particle_name == "gamma") {
				for (G4int i=0;i<2;i++){
					ListOfPrimaryFwdParticles.push_back(theParticleTable->FindParticle(fwd_particle_name));
					ListOfPrimaryAdjParticles.push_back(theParticleTable->FindParticle(adj_particle_name));
				}
			}
			*/
		}
		else {
			if (fwd_ion ){
				ion_name=fwd_ion->GetParticleName();
				G4String adj_ion_name=G4String("adj_") +ion_name;
				ListOfPrimaryFwdParticles.push_back(fwd_ion);
				ListOfPrimaryAdjParticles.push_back(adj_ion);
			}
			else {
				ListOfPrimaryFwdParticles.push_back(0);
				ListOfPrimaryAdjParticles.push_back(0);
				
			}	
		}
	}	
   }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetPrimaryIon(G4ParticleDefinition* adjointIon, G4ParticleDefinition* fwdIon)
{
  fwd_ion = fwdIon;
  adj_ion = adjointIon;
  UpdateListOfPrimaryParticles();
}

