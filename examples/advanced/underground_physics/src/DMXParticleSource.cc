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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
//
// ParticleSource program
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "DMXParticleSource.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"


DMXParticleSource::DMXParticleSource() {

  NumberOfParticlesToBeGenerated = 1;
  particle_definition = NULL;
  G4ThreeVector zero(0., 0., 0.);
  particle_momentum_direction = G4ParticleMomentum(1., 0., 0.);
  particle_energy = 1.0*MeV;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
  particle_charge = 0.0;

  SourcePosType = "Volume";
  Shape = "NULL";
  halfz = 0.;
  Radius = 0.;
  CentreCoords = zero;
  Confine = false;
  VolName = "NULL";

  AngDistType = "iso"; 
  MinTheta = 0.;
  MaxTheta = pi;
  MinPhi = 0.;
  MaxPhi = twopi;

  EnergyDisType = "Mono";
  MonoEnergy = 1*MeV;

  verbosityLevel = 0;

  theMessenger = new DMXParticleSourceMessenger(this);
  gNavigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
}

DMXParticleSource::~DMXParticleSource()
{
  delete theMessenger;
}

void DMXParticleSource::SetPosDisType(G4String PosType) 
{
  SourcePosType = PosType;
}

void DMXParticleSource::SetPosDisShape(G4String shapeType)
{
  Shape = shapeType;
}

void DMXParticleSource::SetCentreCoords(G4ThreeVector coordsOfCentre)
{
  CentreCoords = coordsOfCentre;
}

void DMXParticleSource::SetHalfZ(G4double zhalf)
{
  halfz = zhalf;
}

void DMXParticleSource::SetRadius(G4double radius)
{
  Radius = radius;
}

void DMXParticleSource::ConfineSourceToVolume(G4String Vname)
{
  VolName = Vname;
  if(verbosityLevel == 2) G4cout << VolName << G4endl;

  // checks if selected volume exists
  G4VPhysicalVolume *tempPV      = NULL;
  G4PhysicalVolumeStore *PVStore = 0;
  G4String theRequiredVolumeName = VolName;
  PVStore = G4PhysicalVolumeStore::GetInstance();
  G4int i = 0;
  G4bool found = false;
  if(verbosityLevel == 2) G4cout << PVStore->size() << G4endl;

  // recasting required since PVStore->size() is actually a signed int...
  while (!found && i<(G4int)PVStore->size())
    {
      tempPV = (*PVStore)[i];
      found  = tempPV->GetName() == theRequiredVolumeName;
      if(verbosityLevel == 2)
	G4cout << i << " " << " " << tempPV->GetName() 
	       << " " << theRequiredVolumeName << " " << found << G4endl;
      if (!found)
	{i++;}
    }

  // found = true then the volume exists else it doesnt.
  if(found == true) {
    if(verbosityLevel >= 1)
      G4cout << "Volume " << VolName << " exists" << G4endl;
    Confine = true;
  }
  else if(VolName=="NULL")
    Confine = false;
  else {
    G4cout << " **** Error: Volume does not exist **** " << G4endl;
    G4cout << " Ignoring confine condition" << G4endl;
    VolName = "NULL";
    Confine = false;
  }

}


void DMXParticleSource::SetAngDistType(G4String atype)
{
  AngDistType = atype;
}


void DMXParticleSource::GeneratePointSource()
{
  // Generates Points given the point source.
  if(SourcePosType == "Point")
    particle_position = CentreCoords;
  else
    if(verbosityLevel >= 1)
      G4cout << "Error SourcePosType is not set to Point" << G4endl;
}


void DMXParticleSource::GeneratePointsInVolume()
{
  G4ThreeVector RandPos;
  G4double x=0., y=0., z=0.;
  
  if(SourcePosType != "Volume" && verbosityLevel >= 1)
    G4cout << "Error SourcePosType not Volume" << G4endl;
  
  if(Shape == "Sphere") {
    x = Radius*2.;
    y = Radius*2.;
    z = Radius*2.;
    while(((x*x)+(y*y)+(z*z)) > (Radius*Radius)) {
      x = G4UniformRand();
      y = G4UniformRand();
      z = G4UniformRand();
      
      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
      z = (z*2.*Radius) - Radius;
    }
  }

  else if(Shape == "Cylinder") {
    x = Radius*2.;
    y = Radius*2.;
    while(((x*x)+(y*y)) > (Radius*Radius)) {
      x = G4UniformRand();
      y = G4UniformRand();
      z = G4UniformRand();
      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
      z = (z*2.*halfz) - halfz;
    }
  }
  
  else
    G4cout << "Error: Volume Shape Does Not Exist" << G4endl;

  RandPos.setX(x);
  RandPos.setY(y);
  RandPos.setZ(z);
  particle_position = CentreCoords + RandPos;

}


G4bool DMXParticleSource::IsSourceConfined()
{

  // Method to check point is within the volume specified
  if(Confine == false)
    G4cout << "Error: Confine is false" << G4endl;
  G4ThreeVector null(0.,0.,0.);
  G4ThreeVector *ptr;
  ptr = &null;

  // Check particle_position is within VolName
  G4VPhysicalVolume *theVolume;
  theVolume=gNavigator->LocateGlobalPointAndSetup(particle_position,ptr,true);
  G4String theVolName = theVolume->GetName();
  if(theVolName == VolName) {
    if(verbosityLevel >= 1)
      G4cout << "Particle is in volume " << VolName << G4endl;
    return(true);
  }
  else
    return(false);
}


void DMXParticleSource::SetParticleMomentumDirection
   (G4ParticleMomentum aDirection) {

  particle_momentum_direction =  aDirection.unit();
}


void DMXParticleSource::GenerateIsotropicFlux()
{

  G4double rndm, rndm2;
  G4double px, py, pz;

  G4double sintheta, sinphi, costheta, cosphi;
  rndm = G4UniformRand();
  costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta) - std::cos(MaxTheta));
  sintheta = std::sqrt(1. - costheta*costheta);
  
  rndm2 = G4UniformRand();
  Phi = MinPhi + (MaxPhi - MinPhi) * rndm2; 
  sinphi = std::sin(Phi);
  cosphi = std::cos(Phi);

  px = -sintheta * cosphi;
  py = -sintheta * sinphi;
  pz = -costheta;

  G4double ResMag = std::sqrt((px*px) + (py*py) + (pz*pz));
  px = px/ResMag;
  py = py/ResMag;
  pz = pz/ResMag;

  particle_momentum_direction.setX(px);
  particle_momentum_direction.setY(py);
  particle_momentum_direction.setZ(pz);

  // particle_momentum_direction now holds unit momentum vector.
  if(verbosityLevel >= 2)
    G4cout << "Generating isotropic vector: " << particle_momentum_direction << G4endl;
}


void DMXParticleSource::SetEnergyDisType(G4String DisType)
{
  EnergyDisType = DisType;
}

void DMXParticleSource::SetMonoEnergy(G4double menergy)
{
  MonoEnergy = menergy;
}

void DMXParticleSource::GenerateMonoEnergetic()
{
  particle_energy = MonoEnergy;
}


void DMXParticleSource::SetVerbosity(int vL)
{
  verbosityLevel = vL;
  G4cout << "Verbosity Set to: " << verbosityLevel << G4endl;
}

void DMXParticleSource::SetParticleDefinition
  (G4ParticleDefinition* aParticleDefinition)
{
  particle_definition = aParticleDefinition;
  particle_charge = particle_definition->GetPDGCharge();
}


void DMXParticleSource::GeneratePrimaryVertex(G4Event *evt)
{

  if(particle_definition==NULL) {
    G4cout << "No particle has been defined!" << G4endl;
    return;
  }
  
  // Position
  G4bool srcconf = false;
  G4int LoopCount = 0;
  
  while(srcconf == false)  {
    if(SourcePosType == "Point")
      GeneratePointSource();
    else if(SourcePosType == "Volume")
      GeneratePointsInVolume();
    else {
      G4cout << "Error: SourcePosType undefined" << G4endl;
      G4cout << "Generating point source" << G4endl;
      GeneratePointSource();
    }
    if(Confine == true) {
      srcconf = IsSourceConfined();
      // if source in confined srcconf = true terminating the loop
      // if source isnt confined srcconf = false and loop continues
    }
    else if(Confine == false)
      srcconf = true; // terminate loop
    
    LoopCount++;
    if(LoopCount == 100000) {
      G4cout << "*************************************" << G4endl;
        G4cout << "LoopCount = 100000" << G4endl;
        G4cout << "Either the source distribution >> confinement" << G4endl;
        G4cout << "or any confining volume may not overlap with" << G4endl;
        G4cout << "the source distribution or any confining volumes" << G4endl;
        G4cout << "may not exist"<< G4endl;
        G4cout << "If you have set confine then this will be ignored" <<G4endl;
        G4cout << "for this event." << G4endl;
        G4cout << "*************************************" << G4endl;
        srcconf = true; //Avoids an infinite loop
      }
  }

  // Angular stuff
  if(AngDistType == "iso")
    GenerateIsotropicFlux();
  else if(AngDistType == "direction")
    SetParticleMomentumDirection(particle_momentum_direction);
  else
    G4cout << "Error: AngDistType has unusual value" << G4endl;
  // Energy stuff
  if(EnergyDisType == "Mono")
    GenerateMonoEnergetic();
  else
    G4cout << "Error: EnergyDisType has unusual value" << G4endl;
  
  // create a new vertex
  G4PrimaryVertex* vertex = 
    new G4PrimaryVertex(particle_position,particle_time);

  if(verbosityLevel >= 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;
  // create new primaries and set them to the vertex
  G4double mass =  particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = std::sqrt(energy*energy-mass*mass);
  G4double px = pmom*particle_momentum_direction.x();
  G4double py = pmom*particle_momentum_direction.y();
  G4double pz = pmom*particle_momentum_direction.z();
  
  if(verbosityLevel >= 1){
    G4cout << "Particle name: " 
	   << particle_definition->GetParticleName() << G4endl; 
    G4cout << "       Energy: "<<particle_energy << G4endl;
    G4cout << "     Position: "<<particle_position<< G4endl; 
    G4cout << "    Direction: "<<particle_momentum_direction << G4endl;
    G4cout << " NumberOfParticlesToBeGenerated: "
	   << NumberOfParticlesToBeGenerated << G4endl;
  }


  for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ ) {
    G4PrimaryParticle* particle =
      new G4PrimaryParticle(particle_definition,px,py,pz);
    particle->SetMass( mass );
    particle->SetCharge( particle_charge );
    particle->SetPolarization(particle_polarization.x(),
			      particle_polarization.y(),
			      particle_polarization.z());
    vertex->SetPrimary( particle );
  }
  evt->AddPrimaryVertex( vertex );
  if(verbosityLevel > 1)
    G4cout << " Primary Vetex generated "<< G4endl;   
}




