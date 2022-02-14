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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRActionInitialization.cc
//   Action initialization class
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRActionInitialization.hh"
#include "GRRunAction.hh"
#include "GRPrimaryGeneratorAction.hh"
#include "G4GenericMessenger.hh"

GRActionInitialization::GRActionInitialization()
{
  generatorMsg = new G4GenericMessenger(this,"/gorad/generator/",
                            "primary particle generator selection");

  auto& useParticleGunCmd = generatorMsg->DeclareProperty("useParticleGun",
                useParticleGun, "use Particle Gun");
  useParticleGunCmd.SetStates(G4State_PreInit);
  useParticleGunCmd.SetToBeBroadcasted(false);
  auto& useParticleSourceCmd = generatorMsg->DeclareProperty("useParticleSource",
                useParticleSource, "use General Particle Source");
  useParticleSourceCmd.SetStates(G4State_PreInit);
  useParticleSourceCmd.SetToBeBroadcasted(false);

  filler = new G4TScoreHistFiller<G4AnalysisManager>;
}

GRActionInitialization::~GRActionInitialization()
{
  delete generatorMsg; 
  delete filler;
}

void GRActionInitialization::BuildForMaster() const
{
  SetUserAction(new GRRunAction);
}

void GRActionInitialization::Build() const
{
  SetUserAction(new GRRunAction);
  if(!useParticleGun && !useParticleSource)
  {
    G4ExceptionDescription ed;
    ed << "Neither Particle Gun nor General Particle Source is selected.\n"
       << "No way to generate primary particle!!!\n"
       << "Use /gorad/generator/useParticleGun or /gorad/generator/useParticleSource command.";
    G4Exception("GRActionInitialization::Build()","GORAD0001",FatalException,ed);
  }
  else
  { SetUserAction(new GRPrimaryGeneratorAction(useParticleGun,useParticleSource)); }
}  

