// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTGunGenerator.cc,v 1.1 2003-11-07 21:30:29 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "NTSTGunGenerator.hh"
#include "NTSTGunMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTGunGenerator::NTSTGunGenerator()
: dMeanVertex(0.,0.,0.), dRmsVertex (0.,0.,0.), dPolarization(0.,0.,0.),
  dPlow(1.),  dPhigh(1.), dT0(0.), dCoslow(-1.), dCoshigh(+1.0),       
  dPhilow(0.),dPhihigh(2*pi), dN(1)
{
  messenger = new NTSTGunMessenger(this);
  MeanVertex = dMeanVertex;
  RmsVertex  = dRmsVertex;
  Polarization = dPolarization;
  Plow  = dPlow;
  Phigh = dPhigh;
  Coslow  = dCoslow;
  Coshigh = dCoshigh;
  Philow  = dPhilow;
  Phihigh = dPhihigh;
  T0 = dT0;
  N = dN;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTGunGenerator::~NTSTGunGenerator()
{
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
NTSTGunGenerator::GeneratePrimaryVertex(G4Event* anEvent)
{
  // this function is called at the begining of event
  // 
  // generate primary vertex
  G4double x = MeanVertex.x();
  if (RmsVertex.x()!=0) x += RmsVertex.x()*Gauss();
  G4double y = MeanVertex.y();
  if (RmsVertex.y()!=0) y += RmsVertex.y()*Gauss();
  G4double z = MeanVertex.z();
  if (RmsVertex.z()!=0) z += RmsVertex.z()*Gauss();
  G4double t=0*ns;

  G4PrimaryVertex* vertex = new G4PrimaryVertex(x, y, z, t);

  for (int ipart=0; ipart<GetNumberOfParticles(); ipart++){
    //
    // generate random direction (modifiy later to select an angular range)
    G4double cth = 2.*(G4UniformRand()-0.5);
    G4double sth = sqrt(1.-cth*cth);
    G4double phi   = 2*pi*G4UniformRand();
    G4double cfi = cos(phi);
    G4double sfi = sin(phi);
    G4double p   = (Plow + G4UniformRand()*(Phigh - Plow))*GeV;  
    //
    G4PrimaryParticle* part = new G4PrimaryParticle(GetParticleDefinition(),
						    p*sth*cfi, p*sth*sfi, p*cth);
    //
    // add to vertex
    vertex->SetPrimary( part );
  }
  //
  // add vertex to event
  anEvent->AddPrimaryVertex( vertex );
}




