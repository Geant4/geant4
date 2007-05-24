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
// $Id: DetectorConstruction.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - based on BeamTestDetectorConstruction
//                                    by T. Aso
//
#include "DetectorConstruction.hh"

#include "ConfigData.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4UnitsTable.hh"
#include "Materials.hh"
#include "RingParameterisation.hh"

// Configuration data
using namespace ConfigData;

DetectorConstruction::DetectorConstruction()
  :fLogicalWorld(0)
  ,fDetector(new G4MultiFunctionalDetector(ConfigData::GetDetectorName())) 
{
  Materials::Initialise();
}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  // Material definitions
  Materials::Initialise();
  
  // Construct world
  fLogicalWorld = new G4LogicalVolume(new G4Box("World", 2.0*m, 2.0*m, 2.0*m), 
				      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
				      "World");                                
  
  G4VPhysicalVolume* world = new G4PVPlacement(0,		
					       G4ThreeVector(), 
					       fLogicalWorld,	
					       "World",		
					       0,		
					       false,		
					       0);		
  
  G4SDManager::GetSDMpointer()->AddNewDetector(fDetector);

  // Configuration output file
  std::ofstream& out = ConfigData::GetConfigFile();

  BuildGeometry(out);
  BuildScoringGeometry(out);
  BuildScorers(out);
  
  return world;
}

void DetectorConstruction::BuildGeometry(std::ofstream& out)
{
  CreateTubs("Target", GetTargetMaterial(), GetTargetDistance(), 0.0, GetTargetRadius(), out);
  CreateTubs("ChamberWindow", GetChamberWindowMaterial(), GetChamberWindowDistance(), GetAirGap1Distance(), GetChamberWindowRadius(), out);
  CreateTubs("AirGap1", GetAirGap1Material(), GetAirGap1Distance(), GetMonitorDistance(), GetAirGap1Radius(), out);
  CreateTubs("Monitor", GetMonitorMaterial(), GetMonitorDistance(), GetAirGap2Distance(), GetAirGap1Radius(), out);
  CreateTubs("AirGap2", GetAirGap2Material(), GetAirGap2Distance(), GetBeamWindowDistance(), GetAirGap2Radius(), out);
  CreateTubs("BeamWindow", GetBeamWindowMaterial(), GetBeamWindowDistance(), GetBeamPipeDistance(), GetBeamWindowRadius(), out);
  CreateTubs("BeamPipe", GetBeamPipeMaterial(), GetBeamPipeDistance(), GetBeamPipeEndDistance(), GetBeamPipeRadius(), out);

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction::CreateTubs(G4String name, G4Material* material, 
				      G4double startZ, G4double endZ, G4double radius, std::ofstream& out)
{
  assert (0 != fLogicalWorld);
  
  G4double height = std::fabs(startZ - endZ);
  G4double z = startZ - height/2.0;
  
  G4VSolid* solid = new G4Tubs(name+"_Solid",
			       0.*cm, 
			       radius,
			       height/2.0, 
			       0.*deg, 360.*deg);
  
  G4LogicalVolume* logical = new G4LogicalVolume(solid, material, name+"_Logical");
  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, z),
		    logical, name+"_Physical", fLogicalWorld, false, 0);

  out << "===================================================================" << G4endl;
  out << name << " Summary  "<< G4endl;
  out << " Height  :       " << height/cm << " cm "<< G4endl;
  out << " Radius  :       " << radius/cm << " cm "<< G4endl;
  out << " Z       :       " << z/cm << " cm "<< G4endl;
  out << material <<G4endl;
}

void DetectorConstruction::BuildScoringGeometry(std::ofstream& out) 
{
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  
  // Mother Sphere
  G4VSolid* sphereSolid = new G4Sphere("Sphere_Solid",                   // Name
				       ConfigData::GetScorerMinRadius(), // Min radius
				       ConfigData::GetScorerMaxRadius(), // Max radius
				       0.*deg,                           // Min phi
				       360.*deg,                         // Delta phi
				       0.*deg,                           // Min theta
				       180.*deg);                        // Delta theta

  G4LogicalVolume* sphereLogical = new G4LogicalVolume(sphereSolid, material,"Sphere_Logical");

  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), sphereLogical, "Sphere_Physical", fLogicalWorld, false, 0);

  // Parameterised scoring rings
  G4VSolid* ringsSolid =  new G4Sphere("ScoringRing_Solid",                       // Name
				       ConfigData::GetScorerMinRadius(),   // Min radius
				       ConfigData::GetScorerMaxRadius(),   // Max radius
				       0.*deg,                             // Min phi
				       360.*deg,                           // Delta phi
				       ConfigData::GetScorerMinTheta(),    // Min theta
				       ConfigData::GetScorerDeltaTheta()); // Delta theta
  
  G4LogicalVolume* ringsLogical =  new G4LogicalVolume(ringsSolid, material, "ScoringRing_Logical");


  RingParameterisation* param = new RingParameterisation(ConfigData::GetScorerMinRadius(), 
							 ConfigData::GetScorerMaxRadius(), 
							 ConfigData::GetScorerDeltaTheta(),
							 ConfigData::GetScorerMinTheta());
  
  G4int numberOfReplicas = static_cast<G4int>((ConfigData::GetScorerMaxTheta() 
					       - ConfigData::GetScorerMinTheta())
					      /ConfigData::GetScorerDeltaTheta());
  
  assert (numberOfReplicas > 0) ;
  
  new G4PVParameterised("scorePhys",      // Name
			ringsLogical,     // Logical volume
			sphereLogical,    // Mother volume
			kZAxis,           // Axis
			numberOfReplicas, // Number of replicas
			param);           // Parameterisation
  
  ringsLogical->SetSensitiveDetector(fDetector);
  
  out << "===================================================================" << G4endl;
  out << "Scoring Ring Geometry Summary  "<< G4endl;
  out <<" Minimum radius:   "<< G4BestUnit(ConfigData::GetScorerMinRadius(), "Length") << G4endl;
  out <<" Maximum radius:   "<< G4BestUnit(ConfigData::GetScorerMaxRadius(), "Length") << G4endl;
  out <<" Minimum theta:    "<< ConfigData::GetScorerMinTheta()/deg<<" deg" << G4endl;
  out <<" Maximum theta:    "<< ConfigData::GetScorerMaxTheta()/deg<<" deg" << G4endl;
  out <<" Delta theta:      "<< ConfigData::GetScorerDeltaTheta()/deg<<" deg" << G4endl;
  out <<" Number of strips: "<< numberOfReplicas << G4endl;
  out << material << G4endl;
}

void DetectorConstruction::BuildScorers(std::ofstream& out)
{
  G4double minEnergy = ConfigData::GetScorerMinEnergy();
  G4double maxEnergy = ConfigData::GetScorerMaxEnergy();
  G4double deltaEnergy = ConfigData::GetScorerDeltaEnergy();

  G4int nBins = (maxEnergy - minEnergy)/deltaEnergy;
  assert (nBins > 0);

  for (G4int i=0; i<nBins; ++i) {

    std::ostringstream stream;
    stream<<i;

    G4String name("Scorer"+stream.str());
    G4double lowerEnergy = minEnergy + (i*deltaEnergy);
    G4double upperEnergy = lowerEnergy + deltaEnergy;
    if(ConfigData::GetVerbose()) G4cout<<"Creating scorer: "<<name<<", with energy range: "<<lowerEnergy<<" "<<upperEnergy<<G4endl;

    G4SDParticleWithEnergyFilter* filter = new G4SDParticleWithEnergyFilter("Energy filter", lowerEnergy, upperEnergy);
    filter->add("gamma");

    G4PSSphereSurfaceCurrent* scorer = new G4PSSphereSurfaceCurrent(name, fCurrent_In);
    scorer->SetFilter(filter);

    fDetector->RegisterPrimitive(scorer);    
  }
  
  out << "===================================================================" << G4endl;
  out << "Scorers Summary   "<< G4endl;
  out <<" Minimum Energy:   "<< G4BestUnit(ConfigData::GetScorerMinEnergy(), "Energy") << G4endl;
  out <<" Maximum Energy:   "<< G4BestUnit(ConfigData::GetScorerMaxEnergy(), "Energy") << G4endl;
  out <<" Delta Energy  :   "<< G4BestUnit(ConfigData::GetScorerDeltaEnergy(), "Energy") << G4endl;
  out <<" Number of Bins:   "<< nBins << G4endl;
}
