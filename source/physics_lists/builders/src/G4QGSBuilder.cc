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
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4HadronicParameters.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"

G4QGSBuilder::G4QGSBuilder(const G4String& aName, G4PreCompoundModel*,
			   G4bool quasiel) 
  : G4VHadronModelBuilder(aName), 
    quasielFlag(quasiel)
{}

G4QGSBuilder::~G4QGSBuilder() 
{}                                     

G4HadronicInteraction* G4QGSBuilder::BuildModel()
{
  G4double theMin = G4HadronicParameters::Instance()->GetMinEnergyTransitionQGS_FTF();
  G4double theMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  G4TheoFSGenerator* theModel = new G4TheoFSGenerator(GetName());
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);

  G4QGSModel< G4QGSParticipants >* theStringModel = 
    new G4QGSModel< G4QGSParticipants >;
  G4ExcitedStringDecay* theStringDecay = 
    new G4ExcitedStringDecay(new G4QGSMFragmentation());
  theStringModel->SetFragmentationModel(theStringDecay);

  theModel->SetHighEnergyGenerator(theStringModel);
  if (quasielFlag)
    {
      theModel->SetQuasiElasticChannel(new G4QuasiElasticChannel());
    } 

  if(GetName() == "QGSB") {
    theModel->SetTransport(new G4BinaryCascade());
  } else {
    theModel->SetTransport(new G4GeneratorPrecompoundInterface());
  }

  return theModel;
}

