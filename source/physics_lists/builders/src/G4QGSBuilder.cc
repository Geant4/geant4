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
// $Id: G4QGSBuilder.cc,v 1.1 2009/10/04 16:29:54 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4QGSBuilder
//
// Author: 28 June 2009 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QGSBuilder.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4ProjectileDiffractiveChannel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4QStringChipsParticleLevelInterface.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4BinaryCascade.hh"


G4QGSBuilder::G4QGSBuilder(const G4String& aName, G4bool quasiel, G4bool diff) 
  : G4VHadronModelBuilder(aName), quasielFlag(quasiel), diffFlag(diff)
{
  theProjectileDiffraction = 0;
  theQuasiElastic = 0;
}

G4QGSBuilder::~G4QGSBuilder() 
{
  delete theProjectileDiffraction;
  delete theQuasiElastic;
  delete theQGStringDecay;
  delete theQGStringModel;
}                                     

G4HadronicInteraction* G4QGSBuilder::BuildModel()
{
  G4TheoFSGenerator* theQGSModel = new G4TheoFSGenerator(GetName());
  theQGStringModel  = new G4QGSModel< G4QGSParticipants >;
  theQGStringDecay  = new G4ExcitedStringDecay(new G4QGSMFragmentation());
  theQGStringModel->SetFragmentationModel(theQGStringDecay);
  theQGSModel->SetHighEnergyGenerator(theQGStringModel);

  if(quasielFlag) {
    theQuasiElastic = new G4QuasiElasticChannel();
    theQGSModel->SetQuasiElasticChannel(theQuasiElastic);
  }
  if ( diffFlag ) {
    theProjectileDiffraction = new G4ProjectileDiffractiveChannel();
    theQGSModel->SetProjectileDiffraction(theProjectileDiffraction);
  } 

  if(GetName() == "QGSP") {
    theQGSModel->SetTransport(new G4GeneratorPrecompoundInterface());
  } else if(GetName() == "QGSC") {
    theQGSModel->SetTransport(new G4QStringChipsParticleLevelInterface());
  } else if(GetName() == "QGSB") {
    theQGSModel->SetTransport(new G4BinaryCascade());
  } else {
    theQGSModel->SetTransport(new G4GeneratorPrecompoundInterface());
  }

  return theQGSModel;
}
