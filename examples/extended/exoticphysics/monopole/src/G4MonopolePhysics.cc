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
// $Id: G4MonopolePhysics.cc,v 1.2 2009/07/15 10:19:47 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4MonopolePhysics
//
// Author:      V.Ivanchenko 13.03.2005
//
// Modified:
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MonopolePhysics.hh"
#include "G4MonopolePhysicsMessenger.hh"

#include "G4Monopole.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4StepLimiter.hh"
#include "G4Transportation.hh"
#include "G4MultipleScattering.hh"
#include "G4mplIonisation.hh"
#include "G4hhIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysics::G4MonopolePhysics(const G4String& nam)
  : G4VPhysicsConstructor(nam)
{
  magCharge = 1;
  elCharge  = 0;
  monopoleMass = 100.*GeV;
  theMessenger = new G4MonopolePhysicsMessenger(this);
}

G4MonopolePhysics::~G4MonopolePhysics()
{
  delete theMessenger;
}

void G4MonopolePhysics::ConstructParticle()
{
  G4Monopole::MonopoleDefinition(monopoleMass, magCharge, elCharge);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4MonopolePhysics::ConstructProcess()
{
  if(verboseLevel > 0) {
    G4cout << "G4MonopolePhysics::ConstructProcess" << G4endl;
  }
  G4Monopole* mpl = G4Monopole::Monopole();
  
  G4ProcessManager* pmanager = new G4ProcessManager(mpl);
  mpl->SetProcessManager(pmanager);
  
  // defined monopole parameters and binning

  G4double emax = 10.*TeV;
  G4double magn = mpl->MagneticCharge();
  G4double emin = mpl->GetPDGMass()/20000.;
  if(emin < keV) emin = keV;

  G4int nbin = G4int(std::log10(emin/eV));
  emin       = std::pow(10.,G4double(nbin))*eV;

  nbin = G4int(std::log10(emax/emin));
  if(nbin < 1) nbin = 1;
  nbin *= 10;
  
  pmanager->AddProcess( new G4Transportation(), -1, 0, 0);
  if(magn != 0.0) {
    G4mplIonisation* mplioni = new G4mplIonisation(magn);
    mplioni->SetDEDXBinning(nbin);
    mplioni->SetMinKinEnergy(emin);
    mplioni->SetMaxKinEnergy(emax);
    pmanager->AddProcess(mplioni, -1, 1, 1);
  }
  if(mpl->GetPDGCharge() != 0.0) {
    G4hhIonisation* hhioni = new G4hhIonisation();
    hhioni->SetDEDXBinning(nbin);
    hhioni->SetMinKinEnergy(emin);
    hhioni->SetMaxKinEnergy(emax);
    pmanager->AddProcess(hhioni,  -1, 2, 2);
  }
  pmanager->AddProcess( new G4StepLimiter(),  -1, -1, 3);

}

void G4MonopolePhysics::SetMagneticCharge(G4int val)
{
  magCharge = val;
}

void G4MonopolePhysics::SetElectricCharge(G4int val)
{
  elCharge = val;
}

void G4MonopolePhysics::SetMonopoleMass(G4double mass)
{
  monopoleMass = mass;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

