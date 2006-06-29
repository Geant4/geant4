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
// $Id: NTSTGunGenerator.cc,v 1.5 2006-06-29 18:26:16 gunter Exp $
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
  dPlow(1.),  dPhigh(1.),  dCoslow(-1.), dCoshigh(+1.0),       
  dPhilow(0.),dPhihigh(2*pi), dT0(0.), dN(1)
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
    G4double sth = std::sqrt(1.-cth*cth);
    G4double phi   = 2*pi*G4UniformRand();
    G4double cfi = std::cos(phi);
    G4double sfi = std::sin(phi);
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




