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
// $Id$
//
/////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointCrossSurfChecker
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////

#include "G4AdjointPrimaryGenerator.hh"
#include "G4PhysicalConstants.hh"
#include "G4Event.hh"
#include "G4SingleParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointPosOnPhysVolGenerator.hh" 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPrimaryGenerator::G4AdjointPrimaryGenerator()
  : radius_spherical_source(0.)
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
  	theG4AdjointPosOnPhysVolGenerator->GenerateAPositionOnTheExtSurfaceOfThePhysicalVolume(pos, direction,costh_to_normal);
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
void G4AdjointPrimaryGenerator::SetAdjointPrimarySourceOnAnExtSurfaceOfAVolume(const G4String& volume_name)
{
  theG4AdjointPosOnPhysVolGenerator->DefinePhysicalVolume1(volume_name);
  type_of_adjoint_source ="ExternalSurfaceOfAVolume";
  theSingleParticleSource->GetPosDist()->SetPosDisType("Point");
  theSingleParticleSource->GetAngDist()->SetAngDistType("planar"); 
}

