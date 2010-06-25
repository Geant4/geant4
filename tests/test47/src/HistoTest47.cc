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
#include "globals.hh"
#include "G4ios.hh"

#include "HistoTest47.hh"
#include "G4VParticleChange.hh"

#include "G4AntiProton.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4Proton.hh"

HistoTest47::HistoTest47(std::string namePart, std::string nameMat, 
			 G4double momentum, std::string nameGen) : jobID(-1), 
								   clusterID(-1) {

  setParticle(namePart);
  setTarget(nameMat);
  setMomentum(momentum);
  setGenerator(nameGen);

  mapParticle[G4Proton::Proton()]         = 0;
  mapParticle[G4Neutron::Neutron()]       = 1;
  mapParticle[G4PionPlus::PionPlus()]     = 2;
  mapParticle[G4PionMinus::PionMinus()]   = 3;
  mapParticle[G4KaonPlus::KaonPlus()]     = 4;
  mapParticle[G4KaonMinus::KaonMinus()]   = 5;
  mapParticle[G4AntiProton::AntiProton()] = 6;
  mapParticle[G4PionZero::PionZero()]     = 7;
  unInitialized = true;
}

HistoTest47::~HistoTest47() {}

void HistoTest47::setParticle(std::string namePart) {

  if      (namePart == "pi+") particle = "piplus";
  else if (namePart == "pi-") particle = "piminus";
  else                        particle = namePart; 
  unInitialized = true;
}

G4int HistoTest47::particleType(G4ParticleDefinition* pd) {
  std::map<G4ParticleDefinition*, G4int>::const_iterator it = mapParticle.find(pd);
  if (it != mapParticle.end()) return it->second;
  else                         return -1;
}
