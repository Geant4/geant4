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
// $Id: G4VMscModel.cc,v 1.1 2008-03-07 19:17:53 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VMscModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.03.2008
//
// Modifications:
//
//
// Class Description:
//
// Abstract interface to msc models

// -------------------------------------------------------------------
//

#include "G4VMscModel.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4TransportationManager.hh"
#include "G4SafetyHelper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VMscModel::G4VMscModel(G4double frange, G4double fdtrl, G4double flambdalimit,
			 G4double fgeom,  G4double fskin, G4bool   fsamplez, 
			 G4MscStepLimitType stepAlg, const G4String& nam):
  G4VEmModel(nam), 
  facrange(frange),
  facgeom(fgeom),
  facsafety(0.25),
  skin(fskin),
  dtrl(fdtrl),
  lambdalimit(flambdalimit),
  steppingAlgorithm(stepAlg),
  samplez(fsamplez),
  latDisplasment(true),
  isInitialized(false)
{
  particle      = 0;
  theManager    = G4LossTableManager::Instance(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VMscModel::~G4VMscModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VMscModel::Initialise(const G4ParticleDefinition* p,
			     const G4DataVector&)
{
  if(isInitialized) return;
  // set values of some data members
  SetParticle(p);

  if (pParticleChange)
   fParticleChange = reinterpret_cast<G4ParticleChangeForMSC*>(pParticleChange);
  else
   fParticleChange = new G4ParticleChangeForMSC();

  safetyHelper = G4TransportationManager::GetTransportationManager()
    ->GetSafetyHelper();
  safetyHelper->InitialiseHelper();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


