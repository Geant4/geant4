#include "G4AdjointPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh" 
#include "G4AdjointSimManager.hh" 
#include "G4AdjointPrimaryGenerator.hh"
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPrimaryGeneratorAction::G4AdjointPrimaryGeneratorAction()
{
  
  
  theAdjointPrimaryGenerator= new G4AdjointPrimaryGenerator();

  
  PrimariesConsideredInAdjointSim[G4String("e-")]=false;
  PrimariesConsideredInAdjointSim[G4String("gamma")]=false;
  PrimariesConsideredInAdjointSim[G4String("proton")]=false;
  PrimariesConsideredInAdjointSim[G4String("ion")]=false;


  ListOfPrimaryFwdParticles.clear();
  ListOfPrimaryAdjParticles.clear();
	
 
  
 
  
  last_generated_part_was_adjoint=false;
  index_particle=100000;
  
  ion_name="not_defined";
  fwd_ion = 0;
  adj_ion = 0;
  

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPrimaryGeneratorAction::~G4AdjointPrimaryGeneratorAction()
{delete theAdjointPrimaryGenerator;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{  
  
   
   
   if ( !last_generated_part_was_adjoint ) {
	 
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
   else { //fwd particle equivalent to the last generated adjoint particle ios generated	
  	  G4PrimaryVertex* aPrimVertex = new G4PrimaryVertex();
	  aPrimVertex->SetPosition(pos.x(),pos.y(),pos.z());
  	  aPrimVertex->SetT0(0.);
	  G4PrimaryParticle* aPrimParticle = new G4PrimaryParticle(ListOfPrimaryFwdParticles[index_particle],
  							   -p.x(),-p.y(),-p.z()); 
	
	  aPrimVertex->SetPrimary(aPrimParticle);
 	  anEvent->AddPrimaryVertex(aPrimVertex); 						   
	  last_generated_part_was_adjoint =false;
	  G4AdjointSimManager::GetInstance()->SetAdjointTrackingMode(false);
   }		   
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEmin(G4double val)
{ Emin=val;
  EminIon=val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEmax(G4double val)
{ Emax=val;
  EmaxIon=val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEminIon(G4double val)
{ EminIon=val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::SetEmaxIon(G4double val)
{ EmaxIon=val;
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
void G4AdjointPrimaryGeneratorAction::SetAdjointPrimarySourceOnAnExternalSurfaceOfAVolume(G4String volume_name)
{ type_of_adjoint_source ="ExternalSurfaceOfAVolume";
  theAdjointPrimaryGenerator->SetAdjointPrimarySourceOnAnExternalSurfaceOfAVolume(volume_name); 
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::ConsiderParticleAsPrimary(G4String particle_name)
{ if (PrimariesConsideredInAdjointSim.find(particle_name) != PrimariesConsideredInAdjointSim.end()){
  	PrimariesConsideredInAdjointSim[particle_name]=true;
  }
  UpdateListOfPrimaryParticles();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::NeglectParticleAsPrimary(G4String particle_name)
{ if (PrimariesConsideredInAdjointSim.find(particle_name) != PrimariesConsideredInAdjointSim.end()){
  	PrimariesConsideredInAdjointSim[particle_name]= false;
  }
  UpdateListOfPrimaryParticles();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGeneratorAction::UpdateListOfPrimaryParticles()
{   G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
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
{ fwd_ion = fwdIon;
  adj_ion = adjointIon;
  UpdateListOfPrimaryParticles();
}
