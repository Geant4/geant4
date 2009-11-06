#include "G4AdjointPrimaryGenerator.hh"
#include "G4Event.hh"
#include "G4SingleParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointPosOnPhysVolGenerator.hh" 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPrimaryGenerator::G4AdjointPrimaryGenerator()
{
  theSingleParticleSource  = new G4SingleParticleSource();
 
  theSingleParticleSource->GetEneDist()->SetEnergyDisType("Pow");
  theSingleParticleSource->GetEneDist()->SetAlpha(-1.);
  theSingleParticleSource->GetPosDist()->SetPosDisType("Point");
  theSingleParticleSource->GetAngDist()->SetAngDistType("planar");
  
  
  theG4AdjointPosOnPhysVolGenerator = G4AdjointPosOnPhysVolGenerator::GetInstance();
  
 

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPrimaryGenerator::~G4AdjointPrimaryGenerator()
{
  delete theSingleParticleSource;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGenerator::GenerateAdjointPrimaryVertex(G4Event* anEvent,G4ParticleDefinition* adj_part,G4double E1,G4double E2)
{  
  
   
   
   if (type_of_adjoint_source == "ExternalSurfaceOfAVolume") {
  	
	//Generate position and direction relative to the external surface of sensitive volume
  	//-------------------------------------------------------------

  	G4double costh_to_normal;
	G4ThreeVector pos,direction;
  	theG4AdjointPosOnPhysVolGenerator->GenerateAPositionOnTheExternalSurfaceOfThePhysicalVolume(pos, direction,costh_to_normal);
  	if (costh_to_normal <1.e-4) costh_to_normal =1.e-4;
  	theSingleParticleSource->GetAngDist()->SetParticleMomentumDirection(-direction);
  	theSingleParticleSource->GetPosDist()->SetCentreCoords(pos);
   		
   }	
  		
   theSingleParticleSource->GetEneDist()->SetEmin(E1); 
   theSingleParticleSource->GetEneDist()->SetEmax(E2); 	
	
   theSingleParticleSource->SetParticleDefinition(adj_part);
   theSingleParticleSource->GeneratePrimaryVertex(anEvent);
  	
 	   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGenerator::SetSphericalAdjointPrimarySource(G4double radius, G4ThreeVector center_pos)
{ 
  radius_spherical_source = radius;
  center_spherical_source = center_pos;
  type_of_adjoint_source ="Spherical"; 
  theSingleParticleSource->GetPosDist()->SetPosDisType("Surface");
  theSingleParticleSource->GetPosDist()->SetPosDisShape("Sphere");
  theSingleParticleSource->GetPosDist()->SetCentreCoords(center_pos);
  theSingleParticleSource->GetPosDist()->SetRadius(radius);
  theSingleParticleSource->GetAngDist()->SetAngDistType("cos");
  theSingleParticleSource->GetAngDist()->SetMaxTheta(pi);
  theSingleParticleSource->GetAngDist()->SetMinTheta(halfpi);
  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPrimaryGenerator::SetAdjointPrimarySourceOnAnExternalSurfaceOfAVolume(G4String volume_name)
{ theG4AdjointPosOnPhysVolGenerator->DefinePhysicalVolume1(volume_name);
  type_of_adjoint_source ="ExternalSurfaceOfAVolume";
  theSingleParticleSource->GetPosDist()->SetPosDisType("Point");
  theSingleParticleSource->GetAngDist()->SetAngDistType("planar"); 
}
